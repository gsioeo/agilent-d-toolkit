"""Command line entry point for inspection, search, extraction and quantification.

``inspect`` needs no concentrations and reports what was acquired, including
which channel responded most, so the quantifier can be chosen and written into
``analytes.json``. ``extract`` needs the transitions but no calibration.
``quantify`` runs only with a complete configuration. Every stage writes into a
new directory and never touches a ``.d`` dataset or an earlier result.
"""
import argparse
import sys
import traceback
from pathlib import Path

from . import VERSION, massql_adapter, pipeline, reader_validation, report, search


def _rt_range(value):
    parts = [item for item in value.replace(',', ' ').split() if item]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError('give two numbers, for example "14.2 14.7"')
    result = (float(parts[0]), float(parts[1]))
    if result[0] > result[1]:
        raise argparse.ArgumentTypeError('retention-time start must not exceed end')
    return result


def _nonnegative(value):
    result = float(value)
    if result < 0:
        raise argparse.ArgumentTypeError('value must be non-negative')
    return result


def _fail(out, command, error):
    """Record a clear failure without overwriting anything that already exists."""
    sys.stderr.write('%s failed: %s\n' % (command, error))
    target = Path(out) if out else None
    if target is not None and target.is_dir():
        report.write_json(target / 'failure.json',
                          {'command': command, 'error': str(error),
                           'traceback': traceback.format_exc()})
    return 1


def cmd_inspect(args):
    try:
        summary = pipeline.run_inspect(source=args.source, out=args.out,
                                       ingest_dir=args.ingest_dir,
                                       rt_range=args.rt_range,
                                       checksums=args.checksums)
    except Exception as error:
        return _fail(args.out, 'inspect', error)
    print('inspect: %d run(s), %d channel row(s) -> %s'
          % (summary['runs'], summary['channels'], summary['out_dir']))
    for channel_id, count in sorted(summary['strongest_channel_counts'].items(),
                                    key=lambda item: -item[1]):
        print('  strongest in %d run(s): %s' % (count, channel_id))
    return 0


def cmd_extract(args):
    try:
        summary = pipeline.run_extract(batch_path=args.batch, analytes_path=args.analytes,
                                       out=args.out, data_root=args.data_root,
                                       ingest_dir=args.ingest_dir,
                                       checksums=args.checksums,
                                       write_tic=not args.no_tic)
    except Exception as error:
        return _fail(args.out, 'extract', error)
    print('extract: %d trace(s) from %d run(s) -> %s'
          % (summary['traces'], summary['runs'], summary['out_dir']))
    return 0


def cmd_quantify(args):
    try:
        summary = pipeline.run_quantify(batch_path=args.batch, analytes_path=args.analytes,
                                        out=args.out, data_root=args.data_root,
                                        ingest_dir=args.ingest_dir,
                                        checksums=args.checksums,
                                        plots=not args.no_plots)
    except Exception as error:
        return _fail(args.out, 'quantify', error)
    print('quantify: %d result row(s), %d model(s) -> %s'
          % (summary['results'], summary['models'], summary['out_dir']))
    for failure in summary['failures']:
        print('  failure: %s' % failure)
    return 1 if summary['failures'] else 0


def cmd_search_mrm(args):
    criteria = {'precursor_mz': args.q1, 'product_mz': args.q3,
                'collision_energy_ev': args.ce, 'polarity': args.polarity,
                'time_segment_id': args.segment, 'scan_method_id': args.method,
                'precursor_tol_da': args.mz_tolerance,
                'product_tol_da': args.mz_tolerance,
                'ce_tol_ev': args.ce_tolerance}
    try:
        summary = search.run_mrm_search(
            source=args.source, out=args.out, criteria=criteria,
            ingest_dir=args.ingest_dir, rt_range=args.rt_range,
            checksums=args.checksums)
    except Exception as error:
        return _fail(args.out, 'search-mrm', error)
    print('search-mrm: %d hit(s) in %d run(s) -> %s'
          % (summary['hits'], summary['runs'], summary['out_dir']))
    if summary['unmatched_runs']:
        print('  unmatched: %s' % ', '.join(summary['unmatched_runs']))
    if summary['ambiguous_runs']:
        print('  ambiguous: %s' % ', '.join(summary['ambiguous_runs']))
    return 0


def cmd_search_massql(args):
    try:
        query = args.query
        if args.query_file:
            query = Path(args.query_file).read_text(encoding='utf-8').strip()
        if not query:
            raise ValueError('empty_query: MassQL query must not be empty')
        summary = massql_adapter.run_massql_search(
            source=args.source, out=args.out, query=query,
            ingest_dir=args.ingest_dir, checksums=args.checksums)
    except Exception as error:
        return _fail(args.out, 'search-massql', error)
    print('search-massql: %d result(s) from %d run(s) -> %s'
          % (summary['results'], summary['runs'], summary['out_dir']))
    return 0


def cmd_validate_reader(args):
    try:
        summary = reader_validation.run_reader_validation(
            source=args.source, out=args.out, ingest_dir=args.ingest_dir,
            mz_tolerance=args.mz_tolerance, rt_tolerance=args.rt_tolerance,
            intensity_rtol=args.intensity_rtol,
            intensity_atol=args.intensity_atol, checksums=args.checksums)
    except Exception as error:
        return _fail(args.out, 'validate-reader', error)
    print('validate-reader: %s, %d scan(s), %d failure(s) -> %s'
          % (summary['status'], summary['scans'], summary['failures'],
             summary['out_dir']))
    return 0 if summary['status'].startswith('pass') else 1


def build_parser():
    parser = argparse.ArgumentParser(
        prog='python3 -m mrm_quant',
        description='MRM transition extraction and same-batch calibration quantification')
    parser.add_argument('--version', action='version', version='mrm_quant %s' % VERSION)
    parser.add_argument('--ingest-dir', default=None,
                        help='directory holding agilent_d.py (default: repository ingest/)')
    sub = parser.add_subparsers(dest='command', required=True)

    inspect = sub.add_parser('inspect', help='list acquired channels and rank their response')
    inspect.add_argument('--source', required=True, help='a .d directory or a folder of them')
    inspect.add_argument('--out', required=True, help='new output directory')
    inspect.add_argument('--rt-range', type=_rt_range, default=None,
                         help='restrict the response ranking to a retention-time window')
    inspect.add_argument('--checksums', action='store_true',
                         help='record SHA-256 of the raw files in provenance.json')
    inspect.set_defaults(func=cmd_inspect)

    mrm = sub.add_parser('search-mrm', help='search declared MRM channels and traces')
    mrm.add_argument('--source', required=True, help='a .d directory or a folder of them')
    mrm.add_argument('--out', required=True, help='new output directory')
    mrm.add_argument('--q1', type=float, default=None, help='precursor m/z')
    mrm.add_argument('--q3', type=float, default=None, help='product m/z')
    mrm.add_argument('--ce', type=float, default=None, help='collision energy in eV')
    mrm.add_argument('--polarity', type=int, default=None, help='recorded polarity code')
    mrm.add_argument('--segment', type=int, default=None, help='time segment ID')
    mrm.add_argument('--method', type=int, default=None, help='scan method ID')
    mrm.add_argument('--rt-range', type=_rt_range, default=None,
                     help='restrict response summaries to a retention-time window')
    mrm.add_argument('--mz-tolerance', type=_nonnegative, default=0.01,
                     help='Q1/Q3 matching tolerance in Da (default: 0.01)')
    mrm.add_argument('--ce-tolerance', type=_nonnegative, default=0.01,
                     help='collision-energy tolerance in eV (default: 0.01)')
    mrm.add_argument('--checksums', action='store_true',
                     help='record SHA-256 of the raw files in provenance.json')
    mrm.set_defaults(func=cmd_search_mrm)

    massql = sub.add_parser('search-massql', help='run MassQL directly on .d scans')
    massql.add_argument('--source', required=True, help='a .d directory or a folder of them')
    massql.add_argument('--out', required=True, help='new output directory')
    query = massql.add_mutually_exclusive_group(required=True)
    query.add_argument('--query', help='MassQL query text')
    query.add_argument('--query-file', help='UTF-8 file containing one MassQL query')
    massql.add_argument('--checksums', action='store_true',
                        help='record SHA-256 of the raw files in provenance.json')
    massql.set_defaults(func=cmd_search_massql)

    validation = sub.add_parser('validate-reader',
                                help='compare every centroid scan with rainbow')
    validation.add_argument('--source', required=True,
                            help='a .d directory or a folder of them')
    validation.add_argument('--out', required=True, help='new output directory')
    validation.add_argument('--mz-tolerance', type=_nonnegative, default=1e-6,
                            help='absolute m/z tolerance in Da (default: 1e-6)')
    validation.add_argument('--rt-tolerance', type=_nonnegative, default=1e-6,
                            help='absolute RT tolerance in minutes (default: 1e-6)')
    validation.add_argument('--intensity-rtol', type=_nonnegative, default=1e-6,
                            help='relative intensity tolerance (default: 1e-6)')
    validation.add_argument('--intensity-atol', type=_nonnegative, default=1e-9,
                            help='absolute intensity tolerance (default: 1e-9)')
    validation.add_argument('--checksums', action='store_true',
                            help='record SHA-256 of raw files in provenance.json')
    validation.set_defaults(func=cmd_validate_reader)

    for name, help_text, function in (
            ('extract', 'write the native transition traces', cmd_extract),
            ('quantify', 'integrate, fit the calibration and back-calculate', cmd_quantify)):
        command = sub.add_parser(name, help=help_text)
        command.add_argument('--batch', required=True, help='batch.csv')
        command.add_argument('--analytes', required=True, help='analytes.json')
        command.add_argument('--out', required=True, help='new output directory')
        command.add_argument('--data-root', default=None,
                             help='base for relative dataset paths (default: the batch file directory)')
        command.add_argument('--no-checksums', dest='checksums', action='store_false',
                             help='skip hashing the raw files')
        command.set_defaults(func=function, checksums=True)
    sub.choices['extract'].add_argument('--no-tic', action='store_true',
                                        help='skip the MRM sum chromatogram')
    sub.choices['quantify'].add_argument('--no-plots', action='store_true',
                                         help='skip the figures')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
