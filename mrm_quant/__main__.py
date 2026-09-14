"""Command line entry point: inspect, extract, quantify.

``inspect`` is available now; ``extract`` and ``quantify`` are wired in later
tasks. Every path and parameter comes from the command line or a configuration
file, so no batch, transition or concentration is baked in here.
"""
import argparse
import csv
import sys
from pathlib import Path

from . import reader_adapter


def _datasets(source, reader):
    source = Path(source)
    if source.is_dir() and source.name.lower().endswith('.d'):
        return [source]
    return [Path(p) for p in reader.find_datasets(str(source))]


def _write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def cmd_inspect(args):
    reader = reader_adapter.load_reader(args.ingest_dir)
    out = Path(args.out)
    if out.exists():
        sys.stderr.write('unsafe_output: %s already exists\n' % out)
        return 2
    datasets = _datasets(args.source, reader)
    if not datasets:
        sys.stderr.write('no_dataset_found: %s\n' % args.source)
        return 2

    runs, channels = [], []
    for path in sorted(datasets):
        with reader_adapter.RunReader(path, ingest_dir=args.ingest_dir) as run:
            meta = run.metadata()
            runs.append(meta)
            for channel in run.channels:
                row = dict(channel)
                row['run_name'] = meta['run_name']
                row['method_fingerprint'] = meta['method_fingerprint']
                channels.append(row)

    _write_csv(out / 'runs.csv', runs,
               ['run_name', 'sample_name', 'acquired_time', 'instrument',
                'method_name', 'method_fingerprint', 'scan_count',
                'mrm_frame_count', 'ms1_frame_count', 'rt_min_first',
                'rt_min_last', 'intensity_unit', 'dataset_path'])
    _write_csv(out / 'channels.csv', channels,
               ['run_name', 'channel_id', 'channel_type', 'time_segment_id',
                'scan_method_id', 'method_scan_type', 'compound_name',
                'precursor_mz', 'product_mz', 'collision_energy_ev', 'polarity',
                'dwell_ms', 'gain', 'is_istd', 'precursor_resolution',
                'product_resolution', 'mz_low', 'mz_high', 'step_size',
                'frame_count', 'method_fingerprint'])
    print('inspected %d run(s); wrote %s' % (len(runs), out))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(prog='python3 -m mrm_quant',
                                     description='MRM transition extraction and '
                                                 'same-batch calibration quantification')
    parser.add_argument('--ingest-dir', default=None,
                        help='directory holding agilent_d.py (default: repository ingest/)')
    sub = parser.add_subparsers(dest='command', required=True)

    inspect = sub.add_parser('inspect', help='list acquired channels of a source tree')
    inspect.add_argument('--source', required=True, help='a .d directory or a folder of them')
    inspect.add_argument('--out', required=True, help='new output directory')
    inspect.set_defaults(func=cmd_inspect)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
