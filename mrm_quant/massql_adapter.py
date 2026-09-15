"""Optional MassQL bridge built directly from canonical scan records."""
import datetime

from . import qc, reader_adapter, report


FRAME_COLUMNS = ['i', 'i_norm', 'i_tic_norm', 'mz', 'scan', 'rt', 'polarity']
MS2_COLUMNS = FRAME_COLUMNS + ['precmz', 'ms1scan', 'charge']
CONTEXT_COLUMNS = ['run_id', 'scan', 'rt', 'mslevel', 'scan_type', 'precmz',
                   'collision_energy_ev', 'polarity', 'time_segment_id',
                   'scan_method_id', 'cycle_number', 'point_count']


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _dependencies():
    try:
        import pandas as pd
        from massql import msql_engine
    except ImportError as error:
        raise RuntimeError('optional_dependency_missing: install massql and pandas '
                           'to use search-massql') from error
    return pd, msql_engine


def _peak_rows(record, *, ms1scan=None):
    intensities = [float(value) for value in record.get('intensity') or ()]
    masses = [float(value) for value in record.get('product_mz') or ()]
    total = sum(intensities)
    maximum = max(intensities) if intensities else 0.0
    for mz, intensity in zip(masses, intensities):
        row = {'i': intensity,
               'i_norm': intensity / maximum if maximum else 0.0,
               'i_tic_norm': intensity / total if total else 0.0,
               'mz': mz, 'scan': record.get('scan_id'),
               'rt': record.get('rt_min'), 'polarity': record.get('polarity')}
        if record.get('ms_level') == 2:
            row.update(precmz=record.get('precursor_mz'), ms1scan=ms1scan, charge=None)
        yield row


def records_to_dataframes(records, pandas_module=None):
    """Map canonical records to the DataFrames accepted by MassQL."""
    pd = pandas_module or _dependencies()[0]
    ms1_rows, ms2_rows = [], []
    last_ms1 = None
    for record in records:
        if record.get('ms_level') == 1:
            last_ms1 = record.get('scan_id')
            ms1_rows.extend(_peak_rows(record))
        elif record.get('ms_level') == 2:
            ms2_rows.extend(_peak_rows(record, ms1scan=last_ms1))
    return (pd.DataFrame(ms1_rows, columns=FRAME_COLUMNS),
            pd.DataFrame(ms2_rows, columns=MS2_COLUMNS))


def scan_context(run_id, records):
    """Preserve vendor acquisition fields MassQL does not model."""
    return [{
        'run_id': run_id, 'scan': record.get('scan_id'), 'rt': record.get('rt_min'),
        'mslevel': record.get('ms_level'), 'scan_type': record.get('scan_type'),
        'precmz': record.get('precursor_mz'),
        'collision_energy_ev': record.get('collision_energy_ev'),
        'polarity': record.get('polarity'),
        'time_segment_id': record.get('time_segment_id'),
        'scan_method_id': record.get('scan_method_id'),
        'cycle_number': record.get('cycle_number'),
        'point_count': record.get('point_count')}
        for record in records]


def run_massql_search(*, source, out, query, ingest_dir=None, checksums=False):
    """Run one MassQL query against one or more ``.d`` datasets."""
    pd, engine = _dependencies()
    started = _now()
    datasets = reader_adapter.find_datasets(source, ingest_dir=ingest_dir)
    if not datasets:
        raise ValueError('no_dataset_found: %s' % source)
    out_dir = qc.validate_output_path(out, raw_datasets=datasets)
    out_dir.mkdir(parents=True)
    result_rows, result_columns, contexts, runs, seen_run_ids = [], [], [], [], set()
    for dataset in datasets:
        with reader_adapter.RunReader(dataset, ingest_dir=ingest_dir) as run:
            meta = run.metadata()
            meta.update(run_id=run.run_name, dataset_path=dataset, layout=run.layout)
            runs.append(meta)
            run_id = run.run_name
            if run_id.casefold() in seen_run_ids:
                raise ValueError('duplicate_run_id: %s' % run_id)
            seen_run_ids.add(run_id.casefold())
            records = run.records()
        ms1_df, ms2_df = records_to_dataframes(records, pd)
        result = engine.process_query(query, str(dataset), ms1_df=ms1_df, ms2_df=ms2_df)
        for column in list(result.columns):
            if column not in result_columns:
                result_columns.append(column)
        for row in result.to_dict(orient='records'):
            result_rows.append(dict(row, run_id=run_id))
        contexts.extend(scan_context(run_id, records))
    report.write_csv(out_dir / 'massql_results.csv', result_rows,
                     ['run_id'] + result_columns)
    report.write_csv(out_dir / 'scan_context.csv', contexts, CONTEXT_COLUMNS)
    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='search-massql', out_dir=out_dir, runs=runs, checksums=checksums,
        parameters={'started_utc': started, 'source': str(source), 'query': query,
                    'checksums': checksums},
        extra={'result_count': len(result_rows),
               'note': ('MassQL receives unmodified stored intensities. scan_context.csv '
                        'preserves acquisition fields outside the MassQL data model.')}))
    return {'out_dir': out_dir, 'runs': len(runs), 'results': len(result_rows)}
