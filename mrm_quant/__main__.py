"""Command line entry point: inspect, extract, quantify.

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

from . import VERSION, pipeline, report


def _rt_range(value):
    parts = [item for item in value.replace(',', ' ').split() if item]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError('give two numbers, for example "14.2 14.7"')
    return (float(parts[0]), float(parts[1]))


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
