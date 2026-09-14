"""Transition channel selection and separation of MS1 and MRM signals.

A transition is identified by what the instrument actually acquired: scan mode,
precursor, product, collision energy, polarity, time segment and scan method.
The mass tolerances are matching tolerances for those recorded values, not
isolation windows, and they are never allowed to merge two acquired channels.
A requested transition that the method never acquired is an error; falling back
to an MS1 extracted-ion chromatogram at the same nominal mass is not permitted.
"""
import math

from . import reader_adapter

#: Scan types that carry a selected precursor.
MRM_SCAN_TYPES = reader_adapter.MRM_SCAN_TYPES

#: Keys a selector must provide, plus the matching tolerances.
SELECTOR_KEYS = ('precursor_mz', 'product_mz', 'collision_energy_ev', 'polarity',
                 'time_segment_id', 'scan_method_id')
TOLERANCE_KEYS = ('precursor_tol_da', 'product_tol_da', 'ce_tol_ev')

INTENSITY_UNIT = 'stored_intensity'


def _require(selector, key):
    if key not in selector or selector[key] is None:
        raise ValueError('invalid_selector: %s is required' % key)
    return selector[key]


def _within(value, target, tolerance):
    return value is not None and abs(float(value) - float(target)) <= tolerance


def _matching_indices(values, target, tolerance):
    return [i for i, value in enumerate(values or ()) if _within(value, target, tolerance)]


def _is_mrm(record, mrm_scan_types):
    return record.get('scan_type') in mrm_scan_types and record.get('ms_level') == 2


def selector_from_channel(channel, *, precursor_tol_da=0.01, product_tol_da=0.01,
                          ce_tol_ev=0.01):
    """Build a selector for a method channel produced by :mod:`reader_adapter`."""
    if channel.get('channel_type') != 'mrm':
        raise ValueError('invalid_selector: channel %s is not an MRM channel'
                         % channel.get('channel_id'))
    return {'precursor_mz': channel['precursor_mz'],
            'product_mz': channel['product_mz'],
            'collision_energy_ev': channel['collision_energy_ev'],
            'polarity': channel['polarity'],
            'time_segment_id': channel['time_segment_id'],
            'scan_method_id': channel['scan_method_id'],
            'precursor_tol_da': precursor_tol_da,
            'product_tol_da': product_tol_da,
            'ce_tol_ev': ce_tol_ev}


def header_filter(selector, *, mrm_scan_types=MRM_SCAN_TYPES):
    """Predicate on a canonical record header, used to skip unwanted frames."""
    precursor = _require(selector, 'precursor_mz')
    energy = _require(selector, 'collision_energy_ev')
    polarity = selector.get('polarity')
    segment = selector.get('time_segment_id')
    method_id = selector.get('scan_method_id')
    precursor_tol = selector.get('precursor_tol_da', 0.01)
    ce_tol = selector.get('ce_tol_ev', 0.01)

    def keep(head):
        if not _is_mrm(head, mrm_scan_types):
            return False
        if not _within(head.get('precursor_mz'), precursor, precursor_tol):
            return False
        if not _within(head.get('collision_energy_ev'), energy, ce_tol):
            return False
        if polarity is not None and head.get('polarity') != polarity:
            return False
        if segment is not None and head.get('time_segment_id') != segment:
            return False
        if method_id is not None and head.get('scan_method_id') != method_id:
            return False
        return True

    return keep


def extract_transition(records, selector, *, mrm_scan_types=MRM_SCAN_TYPES):
    """Native trace of one acquired transition.

    Returns equal-length lists ``scan_id``, ``rt_min``, ``intensity`` and
    ``point_status``. A frame that does not belong to the channel is absent from
    the output rather than padded with zero; a declared product that is missing
    from a frame is ``None`` with status ``missing_product``; a measured zero is
    kept as zero with status ``ok``.
    """
    product = _require(selector, 'product_mz')
    product_tol = selector.get('product_tol_da', 0.01)
    keep = header_filter(selector, mrm_scan_types=mrm_scan_types)
    selected = [record for record in records if keep(record)]

    declared_seen = False
    for record in selected:
        declared = record.get('declared_product_mz') or []
        hits = _matching_indices(declared, product, product_tol)
        if len(hits) > 1:
            raise ValueError(
                'ambiguous_product: %d declared products of scan %s are within '
                '%g Da of %g' % (len(hits), record.get('scan_id'), product_tol, product))
        if hits:
            declared_seen = True
    if not declared_seen:
        raise ValueError('transition_not_acquired: no acquired channel declares '
                         'product %g for precursor %g' % (product, selector['precursor_mz']))

    scan_id, rt_min, intensity, point_status, cycle = [], [], [], [], []
    previous_rt = None
    for record in selected:
        current_rt = record.get('rt_min')
        if current_rt is None or not math.isfinite(current_rt):
            raise ValueError('nonfinite: scan %s has no usable retention time'
                             % record.get('scan_id'))
        if previous_rt is not None and current_rt <= previous_rt:
            raise ValueError('non_increasing_time: scan %s repeats or precedes '
                             'retention time %r' % (record.get('scan_id'), previous_rt))
        previous_rt = current_rt

        hits = _matching_indices(record.get('product_mz'), product, product_tol)
        if len(hits) > 1:
            raise ValueError(
                'ambiguous_product: %d measured products of scan %s are within '
                '%g Da of %g' % (len(hits), record.get('scan_id'), product_tol, product))
        if hits:
            value = float(record['intensity'][hits[0]])
            status = 'ok' if math.isfinite(value) else 'nonfinite'
            if status != 'ok':
                value = None
        else:
            value, status = None, 'missing_product'

        scan_id.append(record.get('scan_id'))
        rt_min.append(current_rt)
        intensity.append(value)
        point_status.append(status)
        cycle.append(record.get('cycle_number'))

    return {'scan_id': scan_id, 'rt_min': rt_min, 'intensity': intensity,
            'point_status': point_status, 'cycle_number': cycle,
            'selector': dict(selector), 'intensity_unit': INTENSITY_UNIT,
            'point_count': len(scan_id)}


def extract_run(path, selector, *, ingest_dir=None, mrm_scan_types=MRM_SCAN_TYPES,
                **reader_options):
    """Read one dataset and extract a single transition from it."""
    with reader_adapter.RunReader(path, ingest_dir=ingest_dir, **reader_options) as run:
        records = run.records(header_filter(selector, mrm_scan_types=mrm_scan_types))
        result = extract_transition(records, selector, mrm_scan_types=mrm_scan_types)
    result['dataset_path'] = str(path)
    return result


def _sum_intensity(record):
    total = 0.
    for value in record.get('intensity') or ():
        value = float(value)
        if not math.isfinite(value):
            raise ValueError('nonfinite: scan %s carries a non-finite intensity'
                             % record.get('scan_id'))
        total += value
    return total


def separate_tic(records, *, mrm_scan_types=MRM_SCAN_TYPES):
    """MS1 total ion current and MRM channel sum, kept apart.

    ``ms1_tic`` is the summed signal of the MS1 frames. ``mrm_sum`` is the sum of
    the transitions measured in one MRM acquisition group; it is not a full-scan
    total ion current. Interleaved frames are never joined into a single trace,
    and two acquisition groups are never summed together.
    """
    groups = {'ms1_tic': {}, 'mrm_sum': {}}
    output = {}
    for kind in ('ms1_tic', 'mrm_sum'):
        output[kind] = {'scan_id': [], 'rt_min': [], 'intensity': [],
                        'intensity_unit': INTENSITY_UNIT}
    output['ms1_tic']['description'] = 'summed signal of the MS1 frames'
    output['mrm_sum']['description'] = ('sum of the transitions measured in one '
                                        'MRM acquisition group, not a full-scan TIC')

    for record in records:
        if _is_mrm(record, mrm_scan_types):
            kind = 'mrm_sum'
        elif record.get('ms_level') == 1:
            kind = 'ms1_tic'
        else:
            continue
        key = (record.get('time_segment_id'), record.get('scan_method_id'))
        groups[kind][key] = groups[kind].get(key, 0) + 1
        if len(groups[kind]) > 1:
            raise ValueError('multiple_acquisition_groups: %s covers %s; '
                             'grouped output is required' % (kind, sorted(groups[kind])))
        output[kind]['scan_id'].append(record.get('scan_id'))
        output[kind]['rt_min'].append(record.get('rt_min'))
        output[kind]['intensity'].append(_sum_intensity(record))

    for kind in ('ms1_tic', 'mrm_sum'):
        keys = sorted(groups[kind])
        output[kind]['group'] = keys[0] if keys else None
        output[kind]['point_count'] = len(output[kind]['scan_id'])
    return output


def separate_run(path, *, ingest_dir=None, mrm_scan_types=MRM_SCAN_TYPES,
                 **reader_options):
    """MS1 TIC and MRM sum for one dataset."""
    with reader_adapter.RunReader(path, ingest_dir=ingest_dir, **reader_options) as run:
        return separate_tic(run.records(), mrm_scan_types=mrm_scan_types)
