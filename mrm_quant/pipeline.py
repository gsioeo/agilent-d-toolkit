"""The three stages behind the command line: inspect, extract, quantify.

Each stage writes into a new directory, reads the raw datasets read-only, and
records what it did. Nothing is decided implicitly: the quantifier transition,
the retention-time window, the integration settings, the calibration model and
the quality-control limits all come from the configuration files. ``inspect``
only reports what was acquired, including which channel responded most, so that
choice can be written into the configuration and reviewed.
"""
import datetime
from pathlib import Path

from . import calibration, config, integration, qc, reader_adapter, report, transitions

CHANNEL_COLUMNS = ['run_id', 'channel_id', 'channel_type', 'compound_name',
                   'precursor_mz', 'product_mz', 'collision_energy_ev', 'polarity',
                   'time_segment_id', 'scan_method_id', 'dwell_ms', 'gain', 'is_istd',
                   'precursor_resolution', 'product_resolution', 'mz_low', 'mz_high',
                   'step_size', 'frame_count', 'method_fingerprint']
RUN_COLUMNS = ['run_id', 'sample_name', 'acquired_time', 'instrument', 'method_name',
               'method_fingerprint', 'scan_count', 'mrm_frame_count', 'ms1_frame_count',
               'rt_min_first', 'rt_min_last', 'intensity_unit', 'dataset_path']
RESPONSE_COLUMNS = ['run_id', 'channel_id', 'compound_name', 'precursor_mz', 'product_mz',
                    'max_intensity', 'apex_rt_min', 'area', 'peak_status', 'rank_in_run']
INTEGRATION_COLUMNS = ['run_id', 'analyte_id', 'batch_id', 'role', 'channel_id', 'role_note',
                       'selection_status', 'candidates', 'peak_status', 'apex_rt_min',
                       'apex_intensity', 'height_above_baseline', 'area', 'area_unit',
                       'bounds_start_min', 'bounds_end_min', 'baseline', 'point_count',
                       'expected_rt_min', 'rt_tolerance_min', 'response', 'response_mode',
                       'internal_standard_id', 'internal_standard_area',
                       'internal_standard_status']
CALIBRATION_POINT_COLUMNS = ['batch_id', 'analyte_id', 'run_id', 'level_id', 'concentration',
                             'concentration_unit', 'concentration_source', 'response',
                             'response_mode', 'internal_standard_id',
                             'internal_standard_area',
                             'weight', 'back_calculated_concentration', 'bias_pct',
                             'residual', 'included', 'exclusion_reason']
RESULT_COLUMNS = ['run_id', 'analyte_id', 'batch_id', 'role', 'status', 'validated',
                  'reported_concentration', 'vial_concentration', 'concentration_unit',
                  'dilution_factor', 'area', 'area_unit', 'apex_rt_min', 'response',
                  'response_mode', 'internal_standard_id', 'internal_standard_area',
                  'internal_standard_status', 'calibration_equation', 'r2', 'r2_weighted',
                  'calibration_low', 'calibration_high', 'reasons']
QC_COLUMNS = ['run_id', 'analyte_id', 'batch_id', 'role', 'status', 'validated',
              'qualifier_summary', 'blank_status', 'blank_area', 'blank_limit_area',
              'independent_qc_status', 'qc_bias_pct', 'reasons']
TRACE_COLUMNS = ['scan_id', 'cycle_number', 'rt_min', 'intensity', 'point_status']
TIC_COLUMNS = ['scan_id', 'rt_min', 'intensity']


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _mrm_frames(head):
    return head.get('scan_type') in transitions.MRM_SCAN_TYPES and head.get('ms_level') == 2


def _prepare_output(out, raw_datasets):
    path = qc.validate_output_path(out, raw_datasets=raw_datasets)
    path.mkdir(parents=True)
    return path


def _datasets(source, reader):
    source = Path(source)
    if source.is_dir() and source.name.lower().endswith('.d'):
        return [source]
    return [Path(item) for item in reader.find_datasets(str(source))]


# -- inspect ---------------------------------------------------------------

def run_inspect(*, source, out, ingest_dir=None, rt_range=None, checksums=False):
    """List the acquired channels and rank their response, without any concentration."""
    started = _now()
    reader = reader_adapter.load_reader(ingest_dir)
    datasets = _datasets(source, reader)
    if not datasets:
        raise ValueError('no_dataset_found: %s' % source)
    out_dir = _prepare_output(out, datasets)

    runs, channels, responses = [], [], []
    for dataset in sorted(datasets):
        with reader_adapter.RunReader(dataset, ingest_dir=ingest_dir) as run:
            meta = run.metadata()
            meta['run_id'] = meta['run_name']
            meta['layout'] = run.layout
            meta['dataset_path'] = dataset
            runs.append(meta)
            mrm_channels = [channel for channel in run.channels
                            if channel['channel_type'] == 'mrm']
            for channel in run.channels:
                row = dict(channel)
                row['run_id'] = meta['run_id']
                row['method_fingerprint'] = meta['method_fingerprint']
                channels.append(row)
            if not mrm_channels:
                continue
            records = run.records(_mrm_frames)
        summary = []
        for channel in mrm_channels:
            selector = transitions.selector_from_channel(channel)
            trace = transitions.extract_transition(records, selector)
            values = [value for value in trace['intensity'] if value is not None]
            peaks = integration.find_peaks(trace['rt_min'], trace['intensity'],
                                           rt_range=rt_range, max_peaks=1)
            best = max(peaks, key=lambda peak: peak['apex_intensity']) if peaks else None
            summary.append({'run_id': meta['run_id'], 'channel_id': channel['channel_id'],
                            'compound_name': channel['compound_name'],
                            'precursor_mz': channel['precursor_mz'],
                            'product_mz': channel['product_mz'],
                            'max_intensity': max(values) if values else None,
                            'apex_rt_min': best['apex_rt_min'] if best else None,
                            'area': best['area'] if best else None,
                            'peak_status': best['status'] if best else 'no_peak'})
        summary.sort(key=lambda item: item['max_intensity'] or -1, reverse=True)
        for rank, item in enumerate(summary, start=1):
            item['rank_in_run'] = rank
        responses.extend(summary)

    report.write_csv(out_dir / 'runs.csv', runs, RUN_COLUMNS)
    report.write_csv(out_dir / 'channels.csv', channels, CHANNEL_COLUMNS)
    report.write_csv(out_dir / 'channel_response.csv', responses, RESPONSE_COLUMNS)
    strongest = {}
    for item in responses:
        if item['rank_in_run'] == 1:
            strongest[item['channel_id']] = strongest.get(item['channel_id'], 0) + 1
    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='inspect', out_dir=out_dir, runs=runs, checksums=checksums,
        parameters={'started_utc': started, 'source': str(source),
                    'rt_range': list(rt_range) if rt_range else None,
                    'checksums': checksums},
        extra={'strongest_channel_counts': strongest,
               'note': ('The strongest channel is reported for review only. The '
                        'quantifier transition must be written into analytes.json; '
                        'quantify never chooses a channel by response.')}))
    return {'out_dir': out_dir, 'runs': len(runs), 'channels': len(channels),
            'strongest_channel_counts': strongest}


# -- shared measurement ----------------------------------------------------

def _open_runs(configuration, ingest_dir):
    """Read each dataset once: metadata, fingerprint and its MRM frames."""
    measured = {}
    for row in configuration['batch']['rows']:
        run_id = row['run_id']
        if run_id in measured:
            continue
        with reader_adapter.RunReader(row['dataset_path'], ingest_dir=ingest_dir) as run:
            meta = run.metadata()
            meta['run_id'] = run_id
            meta['dataset_path'] = row['dataset_path']
            meta['layout'] = run.layout
            measured[run_id] = {'meta': meta, 'channels': run.channels,
                                'records': run.records(_mrm_frames)}
    fingerprints = {}
    for run_id, item in measured.items():
        batch_ids = {row['batch_id'] for row in configuration['batch']['rows']
                     if row['run_id'] == run_id}
        for batch_id in batch_ids:
            fingerprints.setdefault(batch_id, {})[run_id] = item['meta']['method_fingerprint']
    for batch_id, mapping in fingerprints.items():
        unique = set(mapping.values())
        if len(unique) > 1:
            raise ValueError('method_mismatch: batch %r mixes method fingerprints %s'
                             % (batch_id, sorted(unique)))
    return measured


def _measure(row, analyte, measured):
    """Select and integrate the quantifier peak, and any configured qualifiers."""
    item = measured[row['run_id']]
    records = item['records']
    settings = analyte['integration']
    window = (analyte['expected_rt_min'] - analyte['search_margin_min'],
              analyte['expected_rt_min'] + analyte['search_margin_min'])
    trace = transitions.extract_transition(records, analyte['quantifier'])
    peaks = integration.find_peaks(
        trace['rt_min'], trace['intensity'], rt_range=window,
        baseline=settings['baseline'], max_gap_min=settings['max_gap_min'],
        min_relative_height=settings['min_relative_height'],
        min_absolute_height=settings['min_absolute_height'],
        boundary_fraction=settings['boundary_fraction'],
        min_separation_min=settings['min_separation_min'],
        min_points=settings['min_points'])
    selection = integration.select_peak(peaks, expected_rt_min=analyte['expected_rt_min'],
                                        tolerance_min=analyte['rt_tolerance_min'])
    peak = selection['peak']
    response, response_mode = None, analyte['response']['mode']
    if peak is not None and peak['status'] == 'ok':
        response = calibration.response_value(peak['area'], mode='external') \
            if response_mode == 'external' else None
    qualifiers = []
    for qualifier in analyte['qualifiers']:
        entry = {'channel_id': qualifier['channel_label'], 'status': 'not_evaluated',
                 'area': None, 'ratio': None}
        try:
            qualifier_trace = transitions.extract_transition(records, qualifier)
        except ValueError as error:
            entry['status'] = str(error).split(':', 1)[0]
            qualifiers.append(entry)
            continue
        if peak is None:
            qualifiers.append(entry)
            continue
        integrated = integration.integrate_peak(
            qualifier_trace['rt_min'], qualifier_trace['intensity'],
            bounds_min=peak['bounds_min'], baseline=settings['baseline'],
            max_gap_min=settings['max_gap_min'], min_points=settings['min_points'])
        entry['area'] = integrated['area']
        reference = qualifier.get('reference_ratio')
        tolerance = qualifier.get('relative_tolerance') \
            or analyte['qc']['qualifier_relative_tolerance']
        if reference is None or tolerance is None:
            entry['status'] = 'not_evaluated'
        else:
            checked = qc.check_qualifier(
                peak['area'], integrated['area'], reference_ratio=reference,
                relative_tolerance=tolerance,
                apex_delta_min=None if integrated['apex_rt_min'] is None
                else integrated['apex_rt_min'] - peak['apex_rt_min'],
                rt_tolerance_min=analyte['qc']['qualifier_rt_tolerance_min'])
            entry.update(checked)
            entry['channel_id'] = qualifier['channel_label']
        qualifiers.append(entry)
    return {'row': row, 'analyte': analyte, 'trace': trace, 'peaks': peaks,
            'selection': selection, 'peak': peak, 'response': response,
            'response_mode': response_mode, 'response_error': None,
            'internal_standard_id': None, 'internal_standard_area': None,
            'internal_standard_status': 'not_used' if response_mode == 'external'
            else 'not_evaluated', 'qualifiers': qualifiers,
            'window': window, 'meta': item['meta']}


def _apply_internal_standard(measurement, analytes, measured):
    """Measure and apply the configured same-injection internal standard."""
    analyte = measurement['analyte']
    standard_id = analyte['response'].get('internal_standard_id')
    measurement['internal_standard_id'] = standard_id
    try:
        standard = analytes[standard_id]
    except (KeyError, TypeError):
        measurement['internal_standard_status'] = 'failed'
        measurement['response_error'] = 'unknown_internal_standard: %r' % standard_id
        return
    try:
        standard_measurement = _measure(measurement['row'], standard, measured)
    except ValueError as error:
        measurement['internal_standard_status'] = 'failed'
        measurement['response_error'] = 'internal_standard_failed: %s' % error
        return
    standard_peak = standard_measurement['peak']
    if standard_peak is None or standard_peak.get('status') != 'ok':
        status = standard_measurement['selection'].get('status')
        if standard_peak is not None:
            status = standard_peak.get('status')
        measurement['internal_standard_status'] = 'failed'
        measurement['response_error'] = 'internal_standard_failed: %s' % (status or 'no_peak')
        return
    measurement['internal_standard_area'] = standard_peak.get('area')
    if measurement['peak'] is None or measurement['peak'].get('status') != 'ok':
        # The IS was valid, but there is no target response to ratio against.
        # Preserve the valid IS diagnostic and let target peak status decide.
        measurement['internal_standard_status'] = 'ok'
        return
    try:
        measurement['response'] = calibration.response_value(
            measurement['peak']['area'],
            mode='internal', is_area=measurement['internal_standard_area'])
    except ValueError as error:
        measurement['internal_standard_status'] = 'failed'
        measurement['response_error'] = str(error)
        return
    measurement['internal_standard_status'] = 'ok'


def _integration_row(measurement):
    row, analyte = measurement['row'], measurement['analyte']
    peak = measurement['peak']
    bounds = peak['bounds_min'] if peak else [None, None]
    return {'run_id': row['run_id'], 'analyte_id': row['analyte_id'],
            'batch_id': row['batch_id'], 'role': row['role'],
            'channel_id': analyte['quantifier']['channel_label'],
            'role_note': 'quantifier',
            'selection_status': measurement['selection']['status'],
            'candidates': measurement['selection']['candidates'],
            'peak_status': peak['status'] if peak else None,
            'apex_rt_min': peak['apex_rt_min'] if peak else None,
            'apex_intensity': peak['apex_intensity'] if peak else None,
            'height_above_baseline': peak.get('height_above_baseline') if peak else None,
            'area': peak['area'] if peak else None,
            'area_unit': integration.AREA_UNIT,
            'bounds_start_min': bounds[0], 'bounds_end_min': bounds[1],
            'baseline': analyte['integration']['baseline'],
            'point_count': peak.get('point_count') if peak else None,
            'expected_rt_min': analyte['expected_rt_min'],
            'rt_tolerance_min': analyte['rt_tolerance_min'],
            'response': measurement['response'],
            'response_mode': measurement['response_mode'],
            'internal_standard_id': measurement['internal_standard_id'],
            'internal_standard_area': measurement['internal_standard_area'],
            'internal_standard_status': measurement['internal_standard_status']}


# -- extract ---------------------------------------------------------------

def run_extract(*, batch_path, analytes_path, out, data_root=None, ingest_dir=None,
                checksums=True, write_tic=True):
    """Write the native transition traces named by the configuration."""
    started = _now()
    configuration = config.load_configuration(batch_path, analytes_path,
                                              data_root=data_root,
                                              require_concentrations=False)
    datasets = [row['dataset_path'] for row in configuration['batch']['rows']]
    out_dir = _prepare_output(out, datasets)
    measured = _open_runs(configuration, ingest_dir)

    channels, written = [], []
    for run_id, item in measured.items():
        for channel in item['channels']:
            row = dict(channel)
            row['run_id'] = run_id
            row['method_fingerprint'] = item['meta']['method_fingerprint']
            channels.append(row)
    for row in configuration['batch']['rows']:
        analyte = configuration['analytes']['analytes'][row['analyte_id']]
        item = measured[row['run_id']]
        for label, selector in [('quantifier', analyte['quantifier'])] + \
                [('qualifier', qualifier) for qualifier in analyte['qualifiers']]:
            trace = transitions.extract_transition(item['records'], selector)
            name = '%s_%s.csv' % (row['analyte_id'], selector['channel_label'])
            path = out_dir / 'traces' / row['run_id'] / name
            report.write_csv(path, [
                {'scan_id': scan, 'cycle_number': cycle, 'rt_min': rt,
                 'intensity': value, 'point_status': status}
                for scan, cycle, rt, value, status in zip(
                    trace['scan_id'], trace['cycle_number'], trace['rt_min'],
                    trace['intensity'], trace['point_status'])], TRACE_COLUMNS)
            written.append({'run_id': row['run_id'], 'analyte_id': row['analyte_id'],
                            'channel_role': label, 'channel_id': selector['channel_label'],
                            'points': trace['point_count'], 'path': str(path)})
        if write_tic:
            separated = transitions.separate_tic(item['records'])
            for kind in ('mrm_sum',):
                path = out_dir / 'tic' / row['run_id'] / ('%s.csv' % kind)
                report.write_csv(path, [
                    {'scan_id': scan, 'rt_min': rt, 'intensity': value}
                    for scan, rt, value in zip(separated[kind]['scan_id'],
                                               separated[kind]['rt_min'],
                                               separated[kind]['intensity'])], TIC_COLUMNS)

    report.write_csv(out_dir / 'channels.csv', channels, CHANNEL_COLUMNS)
    report.write_csv(out_dir / 'traces.csv', written,
                     ['run_id', 'analyte_id', 'channel_role', 'channel_id', 'points', 'path'])
    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='extract', out_dir=out_dir, configuration=configuration,
        runs=[item['meta'] for item in measured.values()], checksums=checksums,
        parameters={'started_utc': started, 'data_root': str(data_root or ''),
                    'write_tic': write_tic, 'checksums': checksums},
        extra={'note': ('MS1 frames are excluded from every transition trace; the '
                        'MRM sum is the sum of the acquired transitions of one '
                        'acquisition group, not a full-scan total ion current.'),
               'internal_standard_note': ('Extract writes configured transitions only; '
                                           'quantify applies any configured internal '
                                           'standard in the same injection.')}))
    return {'out_dir': out_dir, 'traces': len(written), 'runs': len(measured)}


# -- quantify --------------------------------------------------------------

def _blank_area(analyte, measurements):
    blank_run = analyte['qc']['blank_run_id']
    for measurement in measurements:
        row = measurement['row']
        if blank_run is not None and row['run_id'] != blank_run:
            continue
        if blank_run is None and row['role'] != 'blank':
            continue
        peak = measurement['peak']
        return {'run_id': row['run_id'], 'area': peak['area'] if peak else None}
    return None


def run_quantify(*, batch_path, analytes_path, out, data_root=None, ingest_dir=None,
                 checksums=True, plots=True):
    """Integrate, fit the same-batch calibration and back-calculate every run."""
    started = _now()
    configuration = config.load_configuration(batch_path, analytes_path,
                                              data_root=data_root)
    datasets = [row['dataset_path'] for row in configuration['batch']['rows']]
    out_dir = _prepare_output(out, datasets)
    measured = _open_runs(configuration, ingest_dir)

    measurements = []
    for row in configuration['batch']['rows']:
        analyte = configuration['analytes']['analytes'][row['analyte_id']]
        measurements.append(_measure(row, analyte, measured))
    for measurement in measurements:
        if measurement['response_mode'] == 'internal':
            _apply_internal_standard(measurement, configuration['analytes']['analytes'],
                                     measured)

    models, model_rows, results, qc_rows, failures = {}, [], [], [], []
    groups = {}
    for measurement in measurements:
        key = (measurement['row']['batch_id'], measurement['row']['analyte_id'])
        groups.setdefault(key, []).append(measurement)

    for key, group in sorted(groups.items()):
        batch_id, analyte_id = key
        analyte = configuration['analytes']['analytes'][analyte_id]
        fit_rows = []
        for measurement in group:
            row = measurement['row']
            if row['role'] != 'calibration':
                continue
            fit_rows.append({'run_id': row['run_id'], 'batch_id': batch_id,
                             'analyte_id': analyte_id, 'role': 'calibration',
                             'level_id': row['level_id'],
                             'method_fingerprint': measurement['meta']['method_fingerprint'],
                             'concentration': row['concentration'],
                             'concentration_unit': row['concentration_unit'],
                             'response': measurement['response'],
                             'response_error': measurement['response_error'],
                             'response_mode': measurement['response_mode'],
                             'internal_standard_id': measurement['internal_standard_id'],
                             'internal_standard_area': measurement['internal_standard_area'],
                             'include': row['include'],
                             'exclusion_reason': row['exclusion_reason'],
                             'concentration_source': row.get('concentration_source')})
        try:
            model = calibration.fit_calibration(
                fit_rows, weighting=analyte['calibration']['weighting'],
                intercept=analyte['calibration']['intercept'])
        except ValueError as error:
            failures.append({'batch_id': batch_id, 'analyte_id': analyte_id,
                             'stage': 'calibration', 'error': str(error)})
            model = None
        models[key] = model
        if model is not None:
            model['response_mode'] = analyte['response']['mode']
            model['internal_standard_id'] = analyte['response'].get('internal_standard_id')
            by_run = {point['run_id']: point for point in model['points']}
            for fit_row in fit_rows:
                point = by_run.get(fit_row['run_id'], {})
                model_rows.append({
                    'batch_id': batch_id, 'analyte_id': analyte_id,
                    'run_id': fit_row['run_id'], 'level_id': fit_row['level_id'],
                    'concentration': fit_row['concentration'],
                    'concentration_unit': fit_row['concentration_unit'],
                    'concentration_source': fit_row['concentration_source'],
                    'response': fit_row['response'], 'weight': point.get('weight'),
                    'response_mode': fit_row['response_mode'],
                    'internal_standard_id': fit_row['internal_standard_id'],
                    'internal_standard_area': fit_row['internal_standard_area'],
                    'back_calculated_concentration': point.get('back_calculated_concentration'),
                    'bias_pct': point.get('bias_pct'), 'residual': point.get('residual'),
                    'included': fit_row['include'],
                    'exclusion_reason': fit_row['exclusion_reason']})

        blank = _blank_area(analyte, group)
        blank_check = qc.check_blank(None if blank is None else blank['area'],
                                    limit_area=analyte['qc']['blank_limit_area'],
                                    blank_run_id=None if blank is None else blank['run_id'])
        for measurement in group:
            row = measurement['row']
            peak = measurement['peak']
            quantification = None
            if model is not None and measurement['response'] is not None:
                try:
                    quantification = calibration.quantify(
                        measurement['response'], model,
                        dilution_factor=row['dilution_factor'], batch_id=batch_id,
                        method_fingerprint=measurement['meta']['method_fingerprint'])
                except ValueError as error:
                    failures.append({'run_id': row['run_id'], 'analyte_id': analyte_id,
                                     'stage': 'quantify', 'error': str(error)})
                    quantification = {'status': 'quantification_failed',
                                      'original_concentration': None,
                                      'vial_concentration': None, 'error': str(error)}
            elif model is not None and measurement['response_error']:
                quantification = {'status': 'internal_standard_failed',
                                  'original_concentration': None,
                                  'vial_concentration': None,
                                  'error': measurement['response_error']}
            independent = None
            qc_bias = None
            if row['role'] == 'qc' and quantification is not None \
                    and quantification.get('vial_concentration') is not None \
                    and row['concentration'] not in (None, 0):
                qc_bias = (quantification['vial_concentration'] / row['concentration']
                           - 1.) * 100.
                tolerance = analyte['qc']['qc_relative_tolerance']
                if tolerance is None:
                    independent = {'status': 'not_evaluated'}
                else:
                    independent = {'status': 'ok' if abs(qc_bias) <= tolerance * 100.
                                   else 'qc_bias_fail'}
            aggregate = qc.aggregate_status(
                peak_selection=measurement['selection'],
                integration=None if peak is None else {'status': peak['status']},
                quantification=quantification,
                qualifiers=[item for item in measurement['qualifiers']
                            if item['status'] not in ('not_evaluated',)],
                blank=blank_check, independent_qc=independent,
                loq_concentration=analyte['qc']['loq_concentration'],
                calibration_failure=None if model is not None else
                next((failure['error'] for failure in failures
                      if failure.get('batch_id') == batch_id
                      and failure.get('analyte_id') == analyte_id
                      and failure.get('stage') == 'calibration'),
                     'model unavailable'))
            results.append({
                'run_id': row['run_id'], 'analyte_id': analyte_id, 'batch_id': batch_id,
                'role': row['role'], 'status': aggregate['status'],
                'validated': aggregate['validated'],
                'reported_concentration': aggregate['reported_concentration'],
                'vial_concentration': None if quantification is None
                else quantification['vial_concentration'],
                'concentration_unit': row['concentration_unit'],
                'dilution_factor': row['dilution_factor'],
                'area': peak['area'] if peak else None,
                'area_unit': integration.AREA_UNIT,
                'apex_rt_min': peak['apex_rt_min'] if peak else None,
                'response': measurement['response'],
                'response_mode': measurement['response_mode'],
                'internal_standard_id': measurement['internal_standard_id'],
                'internal_standard_area': measurement['internal_standard_area'],
                'internal_standard_status': measurement['internal_standard_status'],
                'calibration_equation': None if model is None else model['equation'],
                'r2': None if model is None else model['r2'],
                'r2_weighted': None if model is None else model['r2_weighted'],
                'calibration_low': None if model is None else model['range'][0],
                'calibration_high': None if model is None else model['range'][1],
                'reasons': '; '.join(aggregate['reasons'])})
            qc_rows.append({
                'run_id': row['run_id'], 'analyte_id': analyte_id, 'batch_id': batch_id,
                'role': row['role'], 'status': aggregate['status'],
                'validated': aggregate['validated'],
                'qualifier_summary': '; '.join(
                    '%s=%s' % (item['channel_id'], item['status'])
                    for item in measurement['qualifiers']) or 'none configured',
                'blank_status': blank_check['status'],
                'blank_area': blank_check['blank_area'],
                'blank_limit_area': blank_check['limit_area'],
                'independent_qc_status': 'not_evaluated' if independent is None
                else independent['status'],
                'qc_bias_pct': qc_bias,
                'reasons': '; '.join(aggregate['reasons'])})

    report.write_csv(out_dir / 'integration.csv',
                     [_integration_row(measurement) for measurement in measurements],
                     INTEGRATION_COLUMNS)
    report.write_csv(out_dir / 'calibration_points.csv', model_rows,
                     CALIBRATION_POINT_COLUMNS)
    report.write_csv(out_dir / 'results.csv', results, RESULT_COLUMNS)
    report.write_csv(out_dir / 'qc.csv', qc_rows, QC_COLUMNS)
    report.write_json(out_dir / 'calibration_models.json', {
        'models': [{'batch_id': key[0], 'analyte_id': key[1],
                    'model': {name: value for name, value in model.items()
                              if name != 'points'},
                    'points': model['points']}
                   for key, model in sorted(models.items()) if model is not None],
        'failures': failures})
    if failures:
        report.write_json(out_dir / 'failures.json', {'failures': failures})

    figures = []
    if plots:
        for key, model in sorted(models.items()):
            if model is None:
                continue
            path = out_dir / 'plots' / ('calibration_%s_%s.png' % key)
            if report.plot_calibration(model, path) is not None:
                figures.append(str(path))
        for measurement in measurements:
            row = measurement['row']
            path = out_dir / 'plots' / 'peaks' / ('%s_%s.png' % (row['run_id'],
                                                                row['analyte_id']))
            title = '%s / %s  %s' % (row['run_id'], row['analyte_id'],
                                     measurement['selection']['status'])
            if report.plot_peak(measurement['trace'], measurement['peak'], path,
                                title=title, window=measurement['window']) is not None:
                figures.append(str(path))

    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='quantify', out_dir=out_dir, configuration=configuration,
        runs=[item['meta'] for item in measured.values()], checksums=checksums,
        parameters={'started_utc': started, 'data_root': str(data_root or ''),
                    'plots': plots, 'checksums': checksums},
        extra={'figures': figures, 'failures': failures,
               'internal_standard_note': ('Internal-standard response uses the '
                                           'configured analyte transition in the same '
                                           'injection. External-standard analytes use '
                                           'no internal standard.'),
               'internal_standard_configuration': {
                   analyte_id: {'mode': analyte['response']['mode'],
                                'internal_standard_id': analyte['response'].get(
                                    'internal_standard_id')}
                   for analyte_id, analyte in
                   configuration['analytes']['analytes'].items()},
               'analyte_settings': {analyte_id: {
                   'expected_rt_min': analyte['expected_rt_min'],
                   'rt_tolerance_min': analyte['rt_tolerance_min'],
                   'search_margin_min': analyte['search_margin_min'],
                   'quantifier': analyte['quantifier'],
                   'qualifiers': analyte['qualifiers'],
                   'integration': analyte['integration'],
                   'calibration': analyte['calibration'],
                   'response': analyte['response'], 'qc': analyte['qc'],
                   'identity_confirmed': analyte['identity_confirmed'],
                   'identity_note': analyte['identity_note']}
                   for analyte_id, analyte
                   in configuration['analytes']['analytes'].items()}}))
    return {'out_dir': out_dir, 'results': len(results),
            'models': sum(model is not None for model in models.values()),
            'failures': failures}
