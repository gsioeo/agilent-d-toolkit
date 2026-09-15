"""Optional scan-level comparison with the independent rainbow reader."""
import datetime

from . import qc, reader_adapter, report


COMPARISON_COLUMNS = [
    'run_id', 'scan_index', 'scan_id', 'status', 'toolkit_rt_min', 'rainbow_rt_min',
    'rt_abs_error_min', 'toolkit_peak_count', 'rainbow_peak_count',
    'mz_max_abs_error_da', 'intensity_max_abs_error', 'detail',
]


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _rainbow():
    try:
        import rainbow
    except ImportError as error:
        raise RuntimeError('optional_dependency_missing: install rainbow-api and numpy '
                           'to use validate-reader') from error
    return rainbow


def _centroid_file(directory):
    for datafile in directory.datafiles:
        if datafile.name.lower() == 'mspeak.bin':
            return datafile
    raise ValueError('rainbow_centroid_missing: MSPeak.bin was not returned')


def _close(a, b, *, atol, rtol=0.0):
    return abs(float(a) - float(b)) <= atol + rtol * abs(float(b))


def _integer_projection(toolkit_values, rainbow_values):
    """Whether rainbow retained exactly the integer part of each intensity."""
    return all(a >= 0 and float(b).is_integer() and int(a) == int(b)
               for a, b in zip(toolkit_values, rainbow_values))


def compare_records(toolkit_records, rainbow_file, *, run_id='', mz_tolerance=1e-6,
                    rt_tolerance=1e-6, intensity_rtol=1e-6,
                    intensity_atol=1e-9):
    """Return one diagnostic row per paired or missing scan."""
    rows = []
    total = max(len(toolkit_records), len(rainbow_file.xlabels))
    for index in range(total):
        if index >= len(toolkit_records):
            rows.append({'run_id': run_id, 'scan_index': index, 'status': 'extra_rainbow',
                         'rainbow_rt_min': float(rainbow_file.xlabels[index]),
                         'detail': 'scan exists only in rainbow'})
            continue
        record = toolkit_records[index]
        base = {'run_id': run_id, 'scan_index': index, 'scan_id': record.get('scan_id'),
                'toolkit_rt_min': record.get('rt_min'),
                'toolkit_peak_count': len(record.get('product_mz') or ())}
        if index >= len(rainbow_file.xlabels):
            base.update(status='missing_rainbow', detail='scan exists only in toolkit')
            rows.append(base)
            continue
        rainbow_rt = float(rainbow_file.xlabels[index])
        rainbow_mz, rainbow_i = rainbow_file.scan(index)
        rainbow_mz = [float(value) for value in rainbow_mz]
        rainbow_i = [float(value) for value in rainbow_i]
        toolkit_mz = [float(value) for value in record.get('product_mz') or ()]
        toolkit_i = [float(value) for value in record.get('intensity') or ()]
        rt_error = abs(float(record['rt_min']) - rainbow_rt)
        base.update(rainbow_rt_min=rainbow_rt, rt_abs_error_min=rt_error,
                    rainbow_peak_count=len(rainbow_mz))
        if len(toolkit_mz) != len(rainbow_mz):
            base.update(status='peak_count_mismatch',
                        detail='%d toolkit peaks, %d rainbow peaks' %
                               (len(toolkit_mz), len(rainbow_mz)))
            rows.append(base)
            continue
        mz_error = max((abs(a - b) for a, b in zip(toolkit_mz, rainbow_mz)), default=0.0)
        intensity_error = max((abs(a - b) for a, b in zip(toolkit_i, rainbow_i)),
                              default=0.0)
        base.update(mz_max_abs_error_da=mz_error,
                    intensity_max_abs_error=intensity_error)
        failures = []
        if rt_error > rt_tolerance:
            failures.append('rt')
        if mz_error > mz_tolerance:
            failures.append('mz')
        intensity_differs = any(
            not _close(a, b, atol=intensity_atol, rtol=intensity_rtol)
            for a, b in zip(toolkit_i, rainbow_i))
        quantized = intensity_differs and _integer_projection(toolkit_i, rainbow_i)
        if intensity_differs and not quantized:
            failures.append('intensity')
        if failures:
            status, detail = 'mismatch', ','.join(failures)
        elif quantized:
            status, detail = 'integer_quantized', 'rainbow retained integer intensity'
        else:
            status, detail = 'ok', ''
        base.update(status=status, detail=detail)
        rows.append(base)
    return rows


def run_reader_validation(*, source, out, ingest_dir=None, mz_tolerance=1e-6,
                          rt_tolerance=1e-6, intensity_rtol=1e-6,
                          intensity_atol=1e-9, checksums=False):
    """Compare every centroid scan from the toolkit and rainbow."""
    rainbow = _rainbow()
    started = _now()
    datasets = reader_adapter.find_datasets(source, ingest_dir=ingest_dir)
    if not datasets:
        raise ValueError('no_dataset_found: %s' % source)
    out_dir = qc.validate_output_path(out, raw_datasets=datasets)
    out_dir.mkdir(parents=True)
    rows, runs, seen_run_ids = [], [], set()
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
        independent = rainbow.read(str(dataset), centroid=True, display_precision=8)
        centroid = _centroid_file(independent)
        rows.extend(compare_records(records, centroid, run_id=run_id,
                                    mz_tolerance=mz_tolerance,
                                    rt_tolerance=rt_tolerance,
                                    intensity_rtol=intensity_rtol,
                                    intensity_atol=intensity_atol))
    accepted = ('ok', 'integer_quantized')
    failures = [row for row in rows if row['status'] not in accepted]
    quantized = sum(row['status'] == 'integer_quantized' for row in rows)
    validation_status = ('fail' if failures else
                         'pass_with_integer_quantization' if quantized else 'pass')
    report.write_csv(out_dir / 'reader_comparison.csv', rows, COMPARISON_COLUMNS)
    report.write_json(out_dir / 'reader_comparison.json', {
        'status': validation_status, 'scan_count': len(rows),
        'failure_count': len(failures),
        'scans_by_status': {status: sum(row['status'] == status for row in rows)
                            for status in sorted({row['status'] for row in rows})},
        'failures_by_status': {status: sum(row['status'] == status for row in failures)
                               for status in sorted({row['status'] for row in failures})},
        'tolerances': {'mz_tolerance_da': mz_tolerance,
                       'rt_tolerance_min': rt_tolerance,
                       'intensity_rtol': intensity_rtol,
                       'intensity_atol': intensity_atol}})
    report.write_json(out_dir / 'provenance.json', report.provenance(
        command='validate-reader', out_dir=out_dir, runs=runs, checksums=checksums,
        parameters={'started_utc': started, 'source': str(source),
                    'mz_tolerance_da': mz_tolerance, 'rt_tolerance_min': rt_tolerance,
                    'intensity_rtol': intensity_rtol,
                    'intensity_atol': intensity_atol, 'checksums': checksums},
        extra={'validation_status': validation_status,
               'failure_count': len(failures),
               'integer_quantized_scan_count': quantized,
               'note': ('This is an independent open-reader comparison, not a '
                        'MassHunter vendor-equivalence claim. integer_quantized means '
                        'rainbow retained the integer part while the toolkit retained '
                        'fractional stored intensity.')}))
    return {'out_dir': out_dir, 'runs': len(runs), 'scans': len(rows),
            'failures': len(failures), 'integer_quantized_scans': quantized,
            'status': validation_status}
