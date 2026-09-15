"""Quality control aggregation and output-path protection."""
import os
import tempfile
import unittest
from pathlib import Path

import context

from mrm_quant import qc


def quantification(status='ok', concentration=0.05, vial=0.05):
    return {'status': status, 'original_concentration': concentration,
            'vial_concentration': vial}


class QualifierTests(unittest.TestCase):
    def test_reference_ratio_must_be_positive(self):
        with self.assertRaisesRegex(ValueError, 'invalid_reference_ratio'):
            qc.check_qualifier(100., 40., reference_ratio=0., relative_tolerance=0.2,
                               apex_delta_min=0., rt_tolerance_min=0.05)

    def test_ratio_and_deviation_are_reported(self):
        result = qc.check_qualifier(200., 60., reference_ratio=0.4, relative_tolerance=0.3,
                                    apex_delta_min=0.0, rt_tolerance_min=0.05)
        self.assertAlmostEqual(result['ratio'], 0.3)
        self.assertAlmostEqual(result['relative_deviation'], 0.25)
        self.assertEqual(result['status'], 'ok')


class BlankTests(unittest.TestCase):
    def test_missing_limit_is_not_a_pass(self):
        self.assertEqual(qc.check_blank(10.)['status'], 'not_evaluated')

    def test_blank_above_the_limit_is_contamination(self):
        self.assertEqual(qc.check_blank(10., limit_area=5.)['status'], 'blank_contamination')
        self.assertEqual(qc.check_blank(4., limit_area=5.)['status'], 'ok')


class AggregationTests(unittest.TestCase):
    def test_clean_result_is_reported_but_unvalidated_without_qc(self):
        result = qc.aggregate_status(peak_selection={'status': 'ok'},
                                     integration={'status': 'ok'},
                                     quantification=quantification())
        self.assertEqual(result['status'], 'ok')
        self.assertAlmostEqual(result['reported_concentration'], 0.05)
        self.assertFalse(result['validated'])

    def test_full_checks_make_the_result_validated(self):
        result = qc.aggregate_status(peak_selection={'status': 'ok'},
                                     integration={'status': 'ok'},
                                     quantification=quantification(),
                                     blank={'status': 'ok'},
                                     independent_qc={'status': 'ok'})
        self.assertTrue(result['validated'])
        self.assertEqual(result['status'], 'ok')

    def test_each_failure_blocks_the_concentration_but_keeps_diagnostics(self):
        cases = [({'peak_selection': {'status': 'ambiguous_peak'}}, 'ambiguous_peak'),
                 ({'integration': {'status': 'non_positive_area'}}, 'non_positive_area'),
                 ({'quantification': quantification('below_calibration_range', None, 0.0001)},
                  'below_calibration_range'),
                 ({'qualifiers': [{'status': 'ion_ratio_fail', 'channel_id': 'q1'}]},
                  'ion_ratio_fail'),
                 ({'blank': {'status': 'blank_contamination', 'blank_run_id': 'BLK'}},
                  'blank_contamination'),
                 ({'independent_qc': {'status': 'qc_bias_fail'}}, 'qc_fail')]
        for overrides, expected in cases:
            with self.subTest(status=expected):
                arguments = {'peak_selection': {'status': 'ok'},
                             'integration': {'status': 'ok'},
                             'quantification': quantification()}
                arguments.update(overrides)
                result = qc.aggregate_status(**arguments)
                self.assertEqual(result['status'], expected)
                self.assertIsNone(result['reported_concentration'])
                self.assertTrue(result['reasons'])

    def test_below_validated_loq_is_distinct_from_below_range(self):
        result = qc.aggregate_status(peak_selection={'status': 'ok'},
                                     integration={'status': 'ok'},
                                     quantification=quantification(concentration=0.001),
                                     loq_concentration=0.01)
        self.assertEqual(result['status'], 'below_validated_loq')
        self.assertIsNone(result['reported_concentration'])
        self.assertAlmostEqual(result['vial_concentration'], 0.05)

    def test_calibration_failure_blocks_result_and_preserves_error(self):
        result = qc.aggregate_status(
            peak_selection={'status': 'ok'}, integration={'status': 'ok'},
            quantification=None, calibration_failure='insufficient_points: 0 usable')
        self.assertEqual(result['status'], 'calibration_failed')
        self.assertIsNone(result['reported_concentration'])
        self.assertIn('insufficient_points: 0 usable', '; '.join(result['reasons']))


class OutputPathTests(unittest.TestCase):
    def test_symlink_into_a_raw_dataset_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / 'run.d'
            raw.mkdir()
            link = root / 'link'
            os.symlink(raw, link)
            with self.assertRaisesRegex(ValueError, 'unsafe_output'):
                qc.validate_output_path(link / 'results', raw_datasets=[raw])

    def test_any_dot_d_ancestor_is_refused_even_if_not_listed(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / 'other.d'
            raw.mkdir()
            with self.assertRaisesRegex(ValueError, 'unsafe_output'):
                qc.validate_output_path(raw / 'deep' / 'results', raw_datasets=[])

    def test_new_directory_is_returned_resolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / 'quant_results' / 'run_001'
            self.assertEqual(qc.validate_output_path(target), target.resolve())

    def test_repository_datasets_cannot_be_used_as_output(self):
        dataset = context.dataset_path('TSJ-0907/LQ/STD_S1.d')
        if not dataset.is_dir():
            self.skipTest('Local raw dataset unavailable')
        with self.assertRaisesRegex(ValueError, 'unsafe_output'):
            qc.validate_output_path(dataset / 'quant', raw_datasets=[dataset])


if __name__ == '__main__':
    unittest.main()
