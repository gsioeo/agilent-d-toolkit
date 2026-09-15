"""Focused pipeline contracts for calibration failure and internal standards."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from mrm_quant import pipeline


def _transition(product):
    return {'precursor_mz': 204., 'product_mz': product,
            'collision_energy_ev': 30., 'polarity': 0,
            'time_segment_id': 1, 'scan_method_id': 1}


class PipelineContractTests(unittest.TestCase):
    def test_failed_calibration_blocks_all_results_and_is_not_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for index, concentration in enumerate((1., 2.), 1):
                run = root / ('cal%d.d' % index); run.mkdir()
                rows.append(['CAL%d' % index, str(run), 'B1', 'A', 'calibration',
                             'S%d' % index, str(concentration)])
            run = root / 'sample.d'; run.mkdir()
            rows.append(['SAMPLE', str(run), 'B1', 'A', 'unknown', '', ''])
            batch = root / 'batch.csv'
            with batch.open('w', newline='', encoding='utf-8') as handle:
                writer = csv.writer(handle)
                writer.writerow(['run_id', 'dataset_path', 'batch_id', 'analyte_id',
                                 'role', 'level_id', 'concentration',
                                 'concentration_unit', 'concentration_basis',
                                 'dilution_factor', 'internal_standard_id',
                                 'internal_standard_concentration', 'include',
                                 'exclusion_reason'])
                for values in rows:
                    writer.writerow(values + ['', '', '1', '', '', 'true', ''])
            analytes = root / 'analytes.json'
            analytes.write_text(json.dumps({'analytes': [{
                'analyte_id': 'A', 'expected_rt_min': 1., 'rt_tolerance_min': .1,
                'quantifier': _transition(93.),
                'calibration': {'concentration_unit': 'ng/mL'}}]}), encoding='utf-8')
            measured = {row[0]: {'meta': {'method_fingerprint': 'M1'}} for row in rows}

            def fake_measure(row, analyte, ignored):
                peak = {'status': 'ok', 'area': 10., 'bounds_min': [0., 1.],
                        'apex_rt_min': 1., 'apex_intensity': 10.,
                        'height_above_baseline': 10., 'point_count': 3}
                return {'row': row, 'analyte': analyte, 'trace': {}, 'peaks': [peak],
                        'selection': {'status': 'ok', 'peak': peak, 'candidates': [peak]},
                        'peak': peak, 'response': 10., 'response_mode': 'external',
                        'response_error': None, 'internal_standard_id': None,
                        'internal_standard_area': None,
                        'internal_standard_status': 'not_used', 'qualifiers': [],
                        'window': [0., 2.], 'meta': ignored[row['run_id']]['meta']}

            with mock.patch.object(pipeline, '_open_runs', return_value=measured), \
                    mock.patch.object(pipeline, '_measure', side_effect=fake_measure):
                summary = pipeline.run_quantify(batch_path=batch, analytes_path=analytes,
                                                out=root / 'out', plots=False,
                                                checksums=False)
            self.assertEqual(summary['models'], 0)
            with (root / 'out' / 'results.csv').open(encoding='utf-8', newline='') as handle:
                results = list(csv.DictReader(handle))
            self.assertTrue(all(row['status'] == 'calibration_failed' for row in results))
            self.assertTrue(all(not row['reported_concentration'] for row in results))
            self.assertTrue(all('insufficient_levels' in row['reasons'] for row in results))
            models = json.loads((root / 'out' / 'calibration_models.json').read_text(
                encoding='utf-8'))
            self.assertEqual(models['models'], [])

    def test_valid_is_is_recorded_when_target_peak_is_absent(self):
        analyte = {'response': {'mode': 'internal', 'internal_standard_id': 'IS'}}
        measurement = {'analyte': analyte, 'row': {'run_id': 'R'},
                       'peak': None, 'response': None,
                       'internal_standard_id': None, 'internal_standard_area': None,
                       'internal_standard_status': 'not_evaluated',
                       'response_error': None}
        standard_peak = {'status': 'ok', 'area': 10.}
        standard = {'peak': standard_peak, 'selection': {'status': 'ok'}}
        with mock.patch.object(pipeline, '_measure', return_value=standard):
            pipeline._apply_internal_standard(measurement, {'IS': {}}, {})
        self.assertEqual(measurement['internal_standard_status'], 'ok')
        self.assertEqual(measurement['internal_standard_area'], 10.)
        self.assertIsNone(measurement['response'])
        self.assertIsNone(measurement['response_error'])

    def test_missing_is_peak_is_an_explicit_failure(self):
        analyte = {'response': {'mode': 'internal', 'internal_standard_id': 'IS'}}
        measurement = {'analyte': analyte, 'row': {'run_id': 'R'},
                       'peak': {'status': 'ok', 'area': 20.}, 'response': None,
                       'internal_standard_id': None, 'internal_standard_area': None,
                       'internal_standard_status': 'not_evaluated',
                       'response_error': None}
        standard = {'peak': None, 'selection': {'status': 'no_peak'}}
        with mock.patch.object(pipeline, '_measure', return_value=standard):
            pipeline._apply_internal_standard(measurement, {'IS': {}}, {})
        self.assertEqual(measurement['internal_standard_status'], 'failed')
        self.assertIn('internal_standard_failed', measurement['response_error'])

    def test_external_measurement_records_no_internal_standard(self):
        measurement = {'response_mode': 'external', 'internal_standard_id': None,
                       'internal_standard_area': None,
                       'internal_standard_status': 'not_used'}
        self.assertEqual(measurement['internal_standard_status'], 'not_used')
        self.assertIsNone(measurement['internal_standard_id'])

    def test_internal_standard_ratio_is_used_and_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            runs = []
            rows = []
            for index, concentration in enumerate((1., 2., 4.), 1):
                run = root / ('cal%d.d' % index); run.mkdir(); runs.append(run)
                rows.append(['CAL%d' % index, str(run), 'B1', 'A', 'calibration',
                             'S%d' % index, str(concentration)])
            run = root / 'sample.d'; run.mkdir(); runs.append(run)
            rows.append(['SAMPLE', str(run), 'B1', 'A', 'unknown', '', ''])
            batch = root / 'batch.csv'
            with batch.open('w', newline='', encoding='utf-8') as handle:
                writer = csv.writer(handle)
                writer.writerow(['run_id', 'dataset_path', 'batch_id', 'analyte_id',
                                 'role', 'level_id', 'concentration',
                                 'concentration_unit', 'concentration_basis',
                                 'dilution_factor', 'internal_standard_id',
                                 'internal_standard_concentration', 'include',
                                 'exclusion_reason'])
                for values in rows:
                    writer.writerow(values + ['', '', '1', '', '', 'true', ''])
            analytes = root / 'analytes.json'
            analytes.write_text(json.dumps({'analytes': [
                {'analyte_id': 'A', 'expected_rt_min': 1., 'rt_tolerance_min': .1,
                 'quantifier': _transition(93.),
                 'calibration': {'concentration_unit': 'ng/mL'},
                 'response': {'mode': 'internal', 'internal_standard_id': 'IS'}},
                {'analyte_id': 'IS', 'expected_rt_min': 1., 'rt_tolerance_min': .1,
                 'quantifier': _transition(81.)}
            ]}), encoding='utf-8')

            measured = {row[0]: {'meta': {'method_fingerprint': 'M1'}} for row in rows}

            def fake_measure(row, analyte, ignored):
                area = 30. if row['role'] == 'unknown' else 20. * row['concentration']
                if analyte['analyte_id'] == 'IS':
                    area = 10.
                peak = {'status': 'ok', 'area': area, 'bounds_min': [0., 1.],
                        'apex_rt_min': 1., 'apex_intensity': area,
                        'height_above_baseline': area, 'point_count': 3}
                return {'row': row, 'analyte': analyte, 'trace': {}, 'peaks': [peak],
                        'selection': {'status': 'ok', 'peak': peak, 'candidates': [peak]},
                        'peak': peak, 'response': area if analyte['response']['mode'] == 'external' else None,
                        'response_mode': analyte['response']['mode'], 'response_error': None,
                        'internal_standard_id': None, 'internal_standard_area': None,
                        'internal_standard_status': 'not_used', 'qualifiers': [],
                        'window': [0., 2.], 'meta': ignored[row['run_id']]['meta']}

            with mock.patch.object(pipeline, '_open_runs', return_value=measured), \
                    mock.patch.object(pipeline, '_measure', side_effect=fake_measure):
                summary = pipeline.run_quantify(batch_path=batch, analytes_path=analytes,
                                                out=root / 'out', plots=False, checksums=False)
            self.assertEqual(summary['models'], 1)
            with (root / 'out' / 'results.csv').open(encoding='utf-8', newline='') as handle:
                result = {row['run_id']: row for row in csv.DictReader(handle)}['SAMPLE']
            self.assertAlmostEqual(float(result['response']), 3.)
            self.assertEqual(result['internal_standard_id'], 'IS')
            self.assertAlmostEqual(float(result['internal_standard_area']), 10.)
            self.assertEqual(result['internal_standard_status'], 'ok')


if __name__ == '__main__':
    unittest.main()
