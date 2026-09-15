"""Configuration: the batch table, the analyte definitions and the level series.

Two UTF-8 files describe a run: ``batch.csv``, one row per injection and analyte,
and ``analytes.json``, one entry per analyte holding the transition, the
retention-time window, the integration and calibration settings and the quality
control limits. Concentrations may be written out level by level, generated from
a series (top concentration, dilution step, number of levels), or both, in which
case the two must agree.
"""
import csv
import hashlib
import json
from pathlib import Path

BATCH_COLUMNS = ('run_id', 'dataset_path', 'batch_id', 'analyte_id', 'role',
                 'level_id', 'concentration', 'concentration_unit',
                 'concentration_basis', 'dilution_factor', 'internal_standard_id',
                 'internal_standard_concentration', 'include', 'exclusion_reason')
REQUIRED_BATCH_COLUMNS = ('run_id', 'dataset_path', 'batch_id', 'analyte_id', 'role')
ROLES = ('calibration', 'blank', 'qc', 'unknown')
TRUE_WORDS = ('1', 'true', 'yes', 'y', 't')
FALSE_WORDS = ('0', 'false', 'no', 'n', 'f')

TRANSITION_KEYS = ('precursor_mz', 'product_mz', 'collision_energy_ev', 'polarity',
                   'time_segment_id', 'scan_method_id')
DEFAULT_TOLERANCES = {'precursor_tol_da': 0.01, 'product_tol_da': 0.01, 'ce_tol_ev': 0.01}
DEFAULT_INTEGRATION = {'baseline': 'linear_endpoints', 'max_gap_min': None,
                       'min_relative_height': 0.02, 'min_absolute_height': None,
                       'boundary_fraction': 0.05, 'min_separation_min': 0.04,
                       'min_points': 3}
DEFAULT_CALIBRATION = {'weighting': 'none', 'intercept': 'free',
                       'concentration_unit': None, 'series': None}
DEFAULT_RESPONSE = {'mode': 'external', 'internal_standard_id': None}
DEFAULT_QC = {'blank_run_id': None, 'blank_limit_area': None, 'loq_concentration': None,
              'qualifier_relative_tolerance': None, 'qualifier_rt_tolerance_min': None,
              'qc_relative_tolerance': None}
DEFAULT_LEVEL_ID_FORMAT = 'S%d'

#: Relative agreement required between a written concentration and the series.
SERIES_RELATIVE_TOLERANCE = 1e-9


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _number(value, field, run_id, *, allow_none=True):
    if value is None or (isinstance(value, str) and not value.strip()):
        if allow_none:
            return None
        raise ValueError('missing_%s: run %r' % (field, run_id))
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError('invalid_%s: run %r has %r' % (field, run_id, value))


def _boolean(value, field, run_id, default=True):
    if value is None or (isinstance(value, str) and not value.strip()):
        return default
    text = str(value).strip().lower()
    if text in TRUE_WORDS:
        return True
    if text in FALSE_WORDS:
        return False
    raise ValueError('invalid_%s: run %r has %r' % (field, run_id, value))


def concentration_series(*, top_concentration, dilution_step, levels,
                         level_id_format=DEFAULT_LEVEL_ID_FORMAT):
    """Level identifier to concentration, halving (or dividing) each step."""
    top = float(top_concentration)
    step = float(dilution_step)
    count = int(levels)
    if not top > 0:
        raise ValueError('invalid_series: top_concentration must be positive, got %r'
                         % top_concentration)
    if not step > 1:
        raise ValueError('invalid_series: dilution_step must be greater than 1, got %r'
                         % dilution_step)
    if count < 1:
        raise ValueError('invalid_series: levels must be at least 1, got %r' % levels)
    series = {}
    for index in range(count):
        try:
            level_id = level_id_format % (index + 1)
        except TypeError:
            raise ValueError('invalid_series: level_id_format %r cannot number a level'
                             % level_id_format)
        series[level_id] = top / step ** index
    return series


def _transition(source, label):
    if not isinstance(source, dict):
        raise ValueError('invalid_transition: %s must be an object' % label)
    transition = {}
    for key in TRANSITION_KEYS:
        if key not in source or source[key] is None:
            raise ValueError('invalid_transition: %s is missing %s' % (label, key))
        transition[key] = source[key]
    for key, default in DEFAULT_TOLERANCES.items():
        transition[key] = source.get(key, default)
    transition['channel_label'] = source.get('channel_label') or (
        '%gto%g' % (float(transition['precursor_mz']), float(transition['product_mz'])))
    return transition


def _merge(defaults, source, label):
    if source is None:
        return dict(defaults)
    if not isinstance(source, dict):
        raise ValueError('invalid_analytes: %s must be an object' % label)
    unknown = set(source) - set(defaults)
    if unknown:
        raise ValueError('invalid_analytes: %s has unknown key(s) %s'
                         % (label, ', '.join(sorted(unknown))))
    merged = dict(defaults)
    merged.update(source)
    return merged


def load_analytes(path):
    """Read and normalise ``analytes.json``."""
    path = Path(path)
    document = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(document, dict) or 'analytes' not in document:
        raise ValueError('invalid_analytes: the file must hold an "analytes" list')
    analytes = {}
    for entry in document['analytes']:
        analyte_id = entry.get('analyte_id')
        if not analyte_id:
            raise ValueError('invalid_analytes: an entry has no analyte_id')
        if analyte_id in analytes:
            raise ValueError('duplicate_analyte: %r appears twice' % analyte_id)
        for key in ('expected_rt_min', 'rt_tolerance_min'):
            if entry.get(key) is None:
                raise ValueError('invalid_analytes: %s is missing %s' % (analyte_id, key))
        calibration = _merge(DEFAULT_CALIBRATION, entry.get('calibration'),
                             '%s.calibration' % analyte_id)
        series = calibration.get('series')
        if series is not None:
            series_levels = concentration_series(
                top_concentration=series.get('top_concentration'),
                dilution_step=series.get('dilution_step'),
                levels=series.get('levels'),
                level_id_format=series.get('level_id_format', DEFAULT_LEVEL_ID_FORMAT))
        else:
            series_levels = None
        qualifiers = []
        for index, qualifier in enumerate(entry.get('qualifiers') or []):
            transition = _transition(qualifier, '%s.qualifiers[%d]' % (analyte_id, index))
            transition['reference_ratio'] = qualifier.get('reference_ratio')
            transition['relative_tolerance'] = qualifier.get('relative_tolerance')
            qualifiers.append(transition)
        analytes[analyte_id] = {
            'analyte_id': analyte_id,
            'identity_confirmed': bool(entry.get('identity_confirmed', False)),
            'identity_note': entry.get('identity_note'),
            'expected_rt_min': float(entry['expected_rt_min']),
            'rt_tolerance_min': float(entry['rt_tolerance_min']),
            'search_margin_min': float(entry.get('search_margin_min',
                                                 3 * float(entry['rt_tolerance_min']))),
            'quantifier': _transition(entry.get('quantifier'), '%s.quantifier' % analyte_id),
            'qualifiers': qualifiers,
            'integration': _merge(DEFAULT_INTEGRATION, entry.get('integration'),
                                  '%s.integration' % analyte_id),
            'calibration': calibration,
            'series_levels': series_levels,
            'response': _merge(DEFAULT_RESPONSE, entry.get('response'),
                               '%s.response' % analyte_id),
            'qc': _merge(DEFAULT_QC, entry.get('qc'), '%s.qc' % analyte_id)}
    if not analytes:
        raise ValueError('invalid_analytes: no analyte is defined')
    for analyte_id, analyte in analytes.items():
        response = analyte['response']
        mode = response.get('mode')
        if mode not in ('external', 'internal'):
            raise ValueError('invalid_response_mode: %r is not external or internal'
                             % mode)
        standard_id = response.get('internal_standard_id')
        if mode == 'external':
            if standard_id is not None:
                raise ValueError('internal_standard_unexpected: %r uses external '
                                 'response but declares %r' % (analyte_id, standard_id))
            continue
        if not standard_id:
            raise ValueError('missing_internal_standard: %r uses internal response '
                             'without internal_standard_id' % analyte_id)
        if standard_id == analyte_id:
            raise ValueError('invalid_internal_standard: %r cannot reference itself'
                             % analyte_id)
        if standard_id not in analytes:
            raise ValueError('unknown_internal_standard: %r refers to %r'
                             % (analyte_id, standard_id))
        if analytes[standard_id]['response'].get('mode') != 'external':
            raise ValueError('invalid_internal_standard: %r must use external '
                             'response' % standard_id)
    return {'analytes': analytes, 'path': str(path), 'sha256': sha256_file(path),
            'version': document.get('version', 'v1'),
            'note': document.get('note')}


def load_batch(path, *, data_root=None, require_datasets=True):
    """Read and validate ``batch.csv``; resolve dataset paths against *data_root*."""
    path = Path(path)
    base = Path(data_root) if data_root else path.parent
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        columns = [name.strip() for name in (reader.fieldnames or [])]
        missing = [name for name in REQUIRED_BATCH_COLUMNS if name not in columns]
        if missing:
            raise ValueError('invalid_batch_columns: missing %s' % ', '.join(missing))
        unknown = [name for name in columns if name not in BATCH_COLUMNS]
        if unknown:
            raise ValueError('invalid_batch_columns: unknown %s' % ', '.join(unknown))
        raw_rows = [dict(row) for row in reader]

    rows, seen, dataset_by_run = [], set(), {}
    for raw in raw_rows:
        run_id = (raw.get('run_id') or '').strip()
        analyte_id = (raw.get('analyte_id') or '').strip()
        if not run_id or not analyte_id:
            raise ValueError('invalid_batch_row: run_id and analyte_id are required')
        key = (run_id, analyte_id)
        if key in seen:
            raise ValueError('duplicate_run_analyte: %r and %r appear twice'
                             % (run_id, analyte_id))
        seen.add(key)
        role = (raw.get('role') or '').strip().lower()
        if role not in ROLES:
            raise ValueError('invalid_role: run %r has %r, expected one of %s'
                             % (run_id, raw.get('role'), ', '.join(ROLES)))
        dataset = (raw.get('dataset_path') or '').strip()
        if not dataset:
            raise ValueError('invalid_batch_row: run %r has no dataset_path' % run_id)
        resolved = Path(dataset)
        if not resolved.is_absolute():
            resolved = base / dataset
        resolved = resolved.expanduser()
        if require_datasets and not resolved.is_dir():
            raise ValueError('dataset_not_found: %s' % resolved)
        previous = dataset_by_run.setdefault(run_id, resolved)
        if previous != resolved:
            raise ValueError('inconsistent_dataset_path: run %r points at %s and %s'
                             % (run_id, previous, resolved))
        include = _boolean(raw.get('include'), 'include', run_id, default=True)
        exclusion_reason = (raw.get('exclusion_reason') or '').strip() or None
        if not include and not exclusion_reason:
            raise ValueError('exclusion_reason_required: run %r is excluded without a reason'
                             % run_id)
        rows.append({
            'run_id': run_id, 'dataset_path': resolved, 'batch_id': (raw.get('batch_id') or '').strip(),
            'analyte_id': analyte_id, 'role': role,
            'level_id': (raw.get('level_id') or '').strip() or None,
            'concentration': _number(raw.get('concentration'), 'concentration', run_id),
            'concentration_unit': (raw.get('concentration_unit') or '').strip() or None,
            'concentration_basis': (raw.get('concentration_basis') or '').strip() or None,
            'dilution_factor': _number(raw.get('dilution_factor'), 'dilution_factor', run_id),
            'internal_standard_id': (raw.get('internal_standard_id') or '').strip() or None,
            'internal_standard_concentration': _number(
                raw.get('internal_standard_concentration'),
                'internal_standard_concentration', run_id),
            'include': include, 'exclusion_reason': exclusion_reason})
    if not rows:
        raise ValueError('invalid_batch_row: the table is empty')
    for row in rows:
        if not row['batch_id']:
            raise ValueError('invalid_batch_row: run %r has no batch_id' % row['run_id'])
        if row['dilution_factor'] is None:
            row['dilution_factor'] = 1.
        elif row['dilution_factor'] <= 0:
            raise ValueError('invalid_dilution: run %r has %r'
                             % (row['run_id'], row['dilution_factor']))
    return {'rows': rows, 'path': str(path), 'sha256': sha256_file(path),
            'data_root': str(base)}


def apply_concentrations(batch, analytes):
    """Fill and cross-check calibration concentrations against the level series."""
    for row in batch['rows']:
        analyte = analytes['analytes'].get(row['analyte_id'])
        if analyte is None:
            raise ValueError('unknown_analyte: run %r refers to %r'
                             % (row['run_id'], row['analyte_id']))
        unit = analyte['calibration'].get('concentration_unit')
        if row['concentration_unit'] is None:
            row['concentration_unit'] = unit
        elif unit is not None and row['concentration_unit'] != unit:
            raise ValueError('unit_mismatch: run %r uses %r, the analyte declares %r'
                             % (row['run_id'], row['concentration_unit'], unit))
        if row['role'] != 'calibration':
            continue
        series = analyte['series_levels']
        expected = None
        if series is not None and row['level_id'] is not None:
            if row['level_id'] not in series:
                raise ValueError('unknown_level: run %r has level %r, the series defines %s'
                                 % (row['run_id'], row['level_id'], ', '.join(sorted(series))))
            expected = series[row['level_id']]
        if row['concentration'] is None:
            if expected is None:
                raise ValueError('missing_concentration: run %r has neither a '
                                 'concentration nor a level in the series' % row['run_id'])
            row['concentration'] = expected
            row['concentration_source'] = 'series'
        else:
            if expected is not None:
                tolerance = SERIES_RELATIVE_TOLERANCE * max(abs(expected), 1e-30)
                if abs(row['concentration'] - expected) > tolerance:
                    raise ValueError('concentration_conflict: run %r is written as %r '
                                     'but the series gives %r'
                                     % (row['run_id'], row['concentration'], expected))
            row['concentration_source'] = 'batch_table'
    return batch


def load_configuration(batch_path, analytes_path, *, data_root=None,
                       require_datasets=True, require_concentrations=True):
    """Load both files, cross-check them and return one configuration object."""
    analytes = load_analytes(analytes_path)
    batch = load_batch(batch_path, data_root=data_root, require_datasets=require_datasets)
    if require_concentrations:
        apply_concentrations(batch, analytes)
    else:
        for row in batch['rows']:
            if row['analyte_id'] not in analytes['analytes']:
                raise ValueError('unknown_analyte: run %r refers to %r'
                                 % (row['run_id'], row['analyte_id']))
    for row in batch['rows']:
        response = analytes['analytes'][row['analyte_id']]['response']
        declared = row.get('internal_standard_id')
        configured = response.get('internal_standard_id')
        if declared is not None and declared != configured:
            raise ValueError('internal_standard_mismatch: run %r declares %r but '
                             'the analyte config declares %r'
                             % (row['run_id'], declared, configured))
    batches = sorted({row['batch_id'] for row in batch['rows']})
    return {'batch': batch, 'analytes': analytes, 'batch_ids': batches,
            'config_hash': hashlib.sha256(
                (batch['sha256'] + analytes['sha256']).encode('ascii')).hexdigest()}
