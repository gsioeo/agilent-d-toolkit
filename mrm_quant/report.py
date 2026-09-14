"""Writers for the result tables, the provenance record and the figures.

CSV keeps full float precision so a value can be read back unchanged; the
figures format for display only. Figures are optional: without matplotlib the
tables are still complete.
"""
import csv
import datetime
import hashlib
import json
import platform
import sys
from pathlib import Path

from . import VERSION, config, reader_adapter

#: Raw files whose checksums identify a dataset for provenance.
CHECKSUM_FILES = ('MSScan.bin', 'MSPeak.bin', 'MSScan.xsd', 'MSTS.xml',
                  'sample_info.xml', 'sequence.xml')

#: Light, plain colours for the figures.
COLOURS = {'point': '#8fb8de', 'line': '#4a7ba7', 'grid': '#dddddd',
           'window': '#eef2f6', 'signal': '#5a5a5a', 'fill': '#cfe0ee'}


def write_csv(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + '\n',
                    encoding='utf-8')
    return path


def dataset_checksums(dataset_path, *, files=CHECKSUM_FILES):
    """SHA-256 of the raw files that define one dataset."""
    acqdata = Path(dataset_path) / 'AcqData'
    digests = {}
    for name in files:
        candidate = acqdata / name
        if candidate.is_file():
            digests[name] = config.sha256_file(candidate)
    method = None
    if acqdata.is_dir():
        for entry in sorted(acqdata.iterdir()):
            if entry.is_dir() and entry.name.lower().endswith('.m'):
                method = entry
                break
    if method is not None:
        for name in ('qqqacqmethod.xml', 'acqmeth.txt'):
            candidate = method / name
            if candidate.is_file():
                digests['%s/%s' % (method.name, name)] = config.sha256_file(candidate)
    return digests


def dependency_versions():
    versions = {'python': sys.version.split()[0], 'platform': platform.platform()}
    for name in ('numpy', 'matplotlib'):
        try:
            module = __import__(name)
        except ImportError:
            versions[name] = None
        else:
            versions[name] = getattr(module, '__version__', 'unknown')
    return versions


def provenance(*, command, parameters, out_dir, configuration=None, runs=(),
               checksums=True, extra=None):
    """Everything needed to repeat the run and to trace every number in it."""
    reader_path = Path(reader_adapter.default_ingest_dir()) / 'agilent_d.py'
    record = {
        'command': command,
        'mrm_quant_version': VERSION,
        'started_utc': parameters.get('started_utc'),
        'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'output_directory': str(out_dir),
        'parameters': parameters,
        'reader': {'path': str(reader_path),
                   'sha256': config.sha256_file(reader_path) if reader_path.is_file() else None},
        'dependencies': dependency_versions(),
        'intensity_unit': 'stored_intensity',
        'area_unit': 'stored_intensity*min',
        'unit_note': ('Intensities are the values stored by the instrument. They are '
                      'not counts or counts per second and the dwell time is not '
                      'applied. No vendor comparison has been performed.'),
        'runs': [],
    }
    if configuration is not None:
        record['configuration'] = {
            'batch_path': configuration['batch']['path'],
            'batch_sha256': configuration['batch']['sha256'],
            'analytes_path': configuration['analytes']['path'],
            'analytes_sha256': configuration['analytes']['sha256'],
            'config_hash': configuration['config_hash'],
            'data_root': configuration['batch']['data_root'],
            'batch_ids': configuration['batch_ids'],
            'exclusions': [{'run_id': row['run_id'], 'analyte_id': row['analyte_id'],
                            'reason': row['exclusion_reason']}
                           for row in configuration['batch']['rows'] if not row['include']]}
    for run in runs:
        entry = {'run_id': run.get('run_id'), 'dataset_path': str(run.get('dataset_path')),
                 'method_fingerprint': run.get('method_fingerprint'),
                 'method_name': run.get('method_name'),
                 'acquired_time': run.get('acquired_time'),
                 'scan_count': run.get('scan_count'),
                 'layout': run.get('layout')}
        if checksums:
            entry['sha256'] = dataset_checksums(run['dataset_path'])
        record['runs'].append(entry)
    if extra:
        record.update(extra)
    return record


def config_digest(text):
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


# -- figures ---------------------------------------------------------------

def _pyplot():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def plot_calibration(model, path, *, title=None):
    """Calibration line with its points, equation and coefficient of determination."""
    plt = _pyplot()
    if plt is None:
        return None
    points = model['points']
    x = [point['concentration'] for point in points]
    y = [point['response'] for point in points]
    figure, axes = plt.subplots(figsize=(6.4, 4.2))
    axes.grid(True, color=COLOURS['grid'], linewidth=0.6)
    axes.plot(x, y, 'o', color=COLOURS['point'], markeredgecolor=COLOURS['line'],
              markersize=5, label='standards')
    low, high = model['range']
    axes.plot([0, high], [model['intercept'], model['slope'] * high + model['intercept']],
              '-', color=COLOURS['line'], linewidth=1.2, label='fit')
    axes.set_xlabel('vial concentration (%s)' % (model.get('concentration_unit') or ''))
    axes.set_ylabel('area (stored_intensity*min)')
    subtitle = '%s   r2 = %s   weighting %s, intercept %s' % (
        model['equation'],
        'n/a' if model['r2'] is None else '%.6f' % model['r2'],
        model['weighting'], model['intercept_mode'])
    axes.set_title('%s\n%s' % (title or '%s / %s' % (model['batch_id'], model['analyte_id']),
                               subtitle), fontsize=9)
    axes.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def plot_peak(trace, peak, path, *, title=None, window=None):
    """One transition around its integration window, with the window shaded."""
    plt = _pyplot()
    if plt is None:
        return None
    rt = trace['rt_min']
    signal = [value if value is not None else float('nan') for value in trace['intensity']]
    figure, axes = plt.subplots(figsize=(6.4, 3.6))
    axes.grid(True, color=COLOURS['grid'], linewidth=0.6)
    if window:
        axes.set_xlim(window)
    axes.plot(rt, signal, '-', color=COLOURS['signal'], linewidth=0.9)
    if peak:
        low, high = peak['bounds_min']
        axes.axvspan(low, high, color=COLOURS['fill'], alpha=0.6)
        axes.plot([peak['apex_rt_min']], [peak['apex_intensity']], 'v',
                  color=COLOURS['line'], markersize=5)
    axes.set_xlabel('retention time (min)')
    axes.set_ylabel('intensity (stored_intensity)')
    if title:
        axes.set_title(title, fontsize=9)
    figure.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path
