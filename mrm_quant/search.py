"""Read-only searches over canonical Agilent scan records."""
import datetime
from pathlib import Path

from . import integration, qc, reader_adapter, report, transitions


MRM_HIT_COLUMNS = [
    'run_id', 'match_status', 'channel_id', 'compound_name', 'precursor_mz',
    'product_mz', 'collision_energy_ev', 'polarity', 'time_segment_id',
    'scan_method_id', 'frame_count', 'trace_point_count', 'max_intensity',
    'apex_rt_min', 'area', 'peak_status', 'trace_path',
]
RUN_COLUMNS = ['run_id', 'sample_name', 'acquired_time', 'instrument', 'method_name',
               'method_fingerprint', 'scan_count', 'mrm_frame_count', 'ms1_frame_count',
               'rt_min_first', 'rt_min_last', 'intensity_unit', 'dataset_path']
TRACE_COLUMNS = ['scan_id', 'cycle_number', 'rt_min', 'intensity', 'point_status']


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _within(value, target, tolerance):
    return target is None or (value is not None and
                              abs(float(value) - float(target)) <= tolerance)


def channel_matches(channel, criteria):
    """Whether one declared MRM channel meets independently optional criteria."""
    if channel.get('channel_type') != 'mrm':
        return False
    if not _within(channel.get('precursor_mz'), criteria.get('precursor_mz'),
                   criteria.get('precursor_tol_da', 0.01)):
        return False
    if not _within(channel.get('product_mz'), criteria.get('product_mz'),
                   criteria.get('product_tol_da', 0.01)):
        return False
    if not _within(channel.get('collision_energy_ev'),
                   criteria.get('collision_energy_ev'),
                   criteria.get('ce_tol_ev', 0.01)):
        return False
    for channel_key, criterion_key in (
            ('polarity', 'polarity'), ('time_segment_id', 'time_segment_id'),
            ('scan_method_id', 'scan_method_id')):
        target = criteria.get(criterion_key)
        if target is not None and channel.get(channel_key) != target:
            return False
    return True


def search_channels(channels, criteria):
    """Return declared MRM channels matching *criteria*, in method order."""
    return [channel for channel in channels if channel_matches(channel, criteria)]


def _trace_rows(trace):
    return [
        {'scan_id': scan, 'cycle_number': cycle, 'rt_min': rt,
         'intensity': value, 'point_status': status}
        for scan, cycle, rt, value, status in zip(
            trace['scan_id'], trace['cycle_number'], trace['rt_min'],
            trace['intensity'], trace['point_status'])
    ]


def run_mrm_search(*, source, out, criteria, ingest_dir=None, rt_range=None,
                   checksums=False):
    """Search declared MRM channels and write their native traces and summaries."""
    started = _now()
    datasets = reader_adapter.find_datasets(source, ingest_dir=ingest_dir)
    if not datasets:
        raise ValueError('no_dataset_found: %s' % source)
    out_dir = qc.validate_output_path(out, raw_datasets=datasets)
    out_dir.mkdir(parents=True)

    hits, runs, unmatched, seen_run_ids = [], [], [], set()
    for dataset in datasets:
        with reader_adapter.RunReader(dataset, ingest_dir=ingest_dir) as run:
            meta = run.metadata()
            meta.update(run_id=run.run_name, dataset_path=dataset, layout=run.layout)
            runs.append(meta)
            run_id = run.run_name
            if run_id.casefold() in seen_run_ids:
                raise ValueError('duplicate_run_id: %s' % run_id)
            seen_run_ids.add(run_id.casefold())
            channels = search_channels(run.channels, criteria)
            if not channels:
                unmatched.append(run_id)
                continue
            records = run.records(lambda head: head.get('scan_type') in
                                  transitions.MRM_SCAN_TYPES and head.get('ms_level') == 2)
        match_status = 'matched' if len(channels) == 1 else 'ambiguous'
        for channel in channels:
            selector = transitions.selector_from_channel(
                channel,
                precursor_tol_da=criteria.get('precursor_tol_da', 0.01),
                product_tol_da=criteria.get('product_tol_da', 0.01),
                ce_tol_ev=criteria.get('ce_tol_ev', 0.01))
            try:
                trace = transitions.extract_transition(records, selector)
            except ValueError as error:
                row = dict(channel)
                row.update(run_id=run_id, match_status=match_status,
                           trace_point_count=0, max_intensity=None,
                           apex_rt_min=None, area=None,
                           peak_status=str(error).split(':', 1)[0], trace_path=None)
                hits.append(row)
                continue
            trace_name = '%s_e%s.csv' % (channel['channel_id'],
                                         channel.get('element_index', 0))
            trace_path = Path('traces') / run_id / trace_name
            report.write_csv(out_dir / trace_path, _trace_rows(trace), TRACE_COLUMNS)
            peaks = integration.find_peaks(trace['rt_min'], trace['intensity'],
                                           rt_range=rt_range, max_peaks=1)
            peak = max(peaks, key=lambda item: item['apex_intensity']) if peaks else None
            values = [value for rt, value in zip(trace['rt_min'], trace['intensity'])
                      if value is not None and
                      (rt_range is None or rt_range[0] <= rt <= rt_range[1])]
            row = dict(channel)
            row.update(run_id=run_id, match_status=match_status,
                       trace_point_count=trace['point_count'],
                       max_intensity=max(values) if values else None,
                       apex_rt_min=peak['apex_rt_min'] if peak else None,
                       area=peak['area'] if peak else None,
                       peak_status=peak['status'] if peak else 'no_peak',
                       trace_path=str(trace_path))
            hits.append(row)

    report.write_csv(out_dir / 'runs.csv', runs, RUN_COLUMNS)
    report.write_csv(out_dir / 'mrm_hits.csv', hits, MRM_HIT_COLUMNS)
    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='search-mrm', out_dir=out_dir, runs=runs, checksums=checksums,
        parameters={'started_utc': started, 'source': str(source),
                    'criteria': criteria,
                    'rt_range': list(rt_range) if rt_range else None,
                    'checksums': checksums},
        extra={'hit_count': len(hits), 'unmatched_runs': unmatched,
               'note': ('Ambiguous matches are reported and never resolved automatically. '
                        'The RT range limits response summaries, not the native trace.')}))
    return {'out_dir': out_dir, 'runs': len(runs), 'hits': len(hits),
            'unmatched_runs': unmatched,
            'ambiguous_runs': sorted({row['run_id'] for row in hits
                                      if row['match_status'] == 'ambiguous'})}
