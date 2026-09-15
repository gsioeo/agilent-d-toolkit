"""Executable specifications for mrm_quant.api.

Default: explicit skips while the package is absent. Set
REQUIRE_MRM_IMPLEMENTATION=1 for the deliberate TDD red gate. Once the package
exists, failures are not hidden. All numbers below are synthetic unless a test
explicitly opens a local .d file.
"""
import copy
import importlib
import importlib.util
import os
import tempfile
import unittest
from pathlib import Path

import context


def transition(**changes):
    result = dict(precursor_mz=204., product_mz=93., collision_energy_ev=30.,
                  polarity=0, time_segment_id=1, scan_method_id=1,
                  precursor_tol_da=0.01, product_tol_da=0.01, ce_tol_ev=0.01)
    result.update(changes)
    return result


def scan(scan_id=1, rt=1., products=(68., 81., 93.), intensities=(1., 2., 10.), **changes):
    result = dict(scan_id=scan_id, rt_min=rt, scan_type=256, ms_level=2,
                  precursor_mz=204., collision_energy_ev=30., polarity=0,
                  time_segment_id=1, scan_method_id=1, cycle_number=scan_id,
                  product_mz=list(products), intensity=list(intensities),
                  declared_product_mz=[68., 81., 93.])
    result.update(changes)
    return result


def standards():
    return [dict(run_id=f'cal{i}', batch_id='B1', analyte_id='A',
                 method_fingerprint='M1', role='calibration', concentration=x,
                 concentration_unit='ng/mL', response=2*x+1,
                 include=True, exclusion_reason=None)
            for i, x in enumerate([1., 2., 4.])]


class FutureContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if importlib.util.find_spec('mrm_quant') is None:
            if os.environ.get('REQUIRE_MRM_IMPLEMENTATION') == '1':
                raise ModuleNotFoundError('Planned mrm_quant.api is intentionally not implemented')
            raise unittest.SkipTest('Planned mrm_quant.api is not implemented; these are pending contracts')
        cls.api = importlib.import_module('mrm_quant.api')

    def test_transition_filters_ms1_other_precursor_ce_polarity_and_segments(self):
        records = [scan(), scan(2, 1.1, ms_level=1, scan_type=1),
                   scan(3, 1.2, precursor_mz=205.), scan(4, 1.3, collision_energy_ev=20.),
                   scan(5, 1.4, polarity=1), scan(6, 1.5, time_segment_id=2),
                   scan(7, 1.6, scan_method_id=2), scan(8, 1.7, scan_type=512)]
        result = self.api.extract_transition(records, transition())
        self.assertEqual(result['scan_id'], [1])
        self.assertEqual(result['rt_min'], [1.])
        self.assertEqual(result['intensity'], [10.])

    def test_product_values_are_matched_by_mz_not_position(self):
        result = self.api.extract_transition([scan(products=(93., 68., 81.), intensities=(10., 1., 2.))], transition())
        self.assertEqual(result['intensity'], [10.])

    def test_declared_but_missing_product_is_missing_not_zero(self):
        result = self.api.extract_transition([scan(products=(68., 81.), intensities=(1., 2.))], transition())
        self.assertEqual(result['intensity'], [None])
        self.assertEqual(result['point_status'], ['missing_product'])

    def test_measured_zero_is_preserved(self):
        result = self.api.extract_transition([scan(intensities=(1., 2., 0.))], transition())
        self.assertEqual(result['intensity'], [0.])
        self.assertEqual(result['point_status'], ['ok'])

    def test_unacquired_transition_cannot_fall_back_to_full_scan(self):
        with self.assertRaisesRegex(ValueError, 'transition_not_acquired'):
            self.api.extract_transition([scan(), scan(2, 1.1, ms_level=1, scan_type=1)], transition(product_mz=189.))

    def test_ambiguous_product_window_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'ambiguous_product'):
            self.api.extract_transition([scan(products=(92.99, 93.01), intensities=(5., 8.),
                                             declared_product_mz=[92.99, 93.01])], transition(product_tol_da=0.03))

    def test_duplicate_times_are_rejected_within_transition(self):
        with self.assertRaisesRegex(ValueError, 'non_increasing_time'):
            self.api.extract_transition([scan(1, 1.), scan(2, 1.)], transition())

    def test_mass_tolerance_does_not_round_to_integer(self):
        result = self.api.extract_transition([scan(products=(93.009,), intensities=(7.,))], transition())
        self.assertEqual(result['intensity'], [7.])
        result = self.api.extract_transition([scan(products=(93.02,), intensities=(7.,))], transition())
        self.assertEqual(result['intensity'], [None])

    def test_optional_real_transition_preserves_only_mrm_times(self):
        path = context.dataset_path('mrm_primary')
        if not path.exists():
            self.skipTest('Local raw fixture unavailable')
        result = self.api.extract_run(path, transition())
        self.assertTrue(result['intensity'])
        self.assertEqual(len(result['scan_id']), len(result['rt_min']))
        self.assertTrue(all(b > a for a, b in zip(result['rt_min'], result['rt_min'][1:])))

    def test_ms1_tic_and_mrm_sum_are_separate(self):
        result = self.api.separate_tic([scan(), scan(2, 1.1, scan_type=1, ms_level=1,
                                                  scan_method_id=2, intensities=(100., 200., 300.))])
        self.assertEqual(result['ms1_tic']['scan_id'], [2])
        self.assertEqual(result['ms1_tic']['intensity'], [600.])
        self.assertEqual(result['mrm_sum']['scan_id'], [1])
        self.assertEqual(result['mrm_sum']['intensity'], [13.])

    def test_irregular_time_linear_baseline_area_has_minute_units(self):
        # Baseline [10,11,14,16], excess [0,2,2,0]; trapezoids: 1+6+2=9.
        result = self.api.integrate_peak([0., 1., 4., 6.], [10., 13., 16., 16.],
                                         bounds_min=(0., 6.), baseline='linear_endpoints', max_gap_min=4.)
        self.assertEqual(result['status'], 'ok')
        self.assertAlmostEqual(result['area'], 9.)
        self.assertEqual(result['area_unit'], 'stored_intensity*min')

    def test_integration_interpolates_window_endpoints(self):
        result = self.api.integrate_peak([0., 1., 2.], [0., 2., 0.],
                                         bounds_min=(0.5, 1.5), baseline='none', max_gap_min=2.)
        self.assertAlmostEqual(result['area'], 1.5)

    def test_integration_rejects_missing_nonfinite_times_and_large_gaps(self):
        for rt, signal, code in [([0., 1., 2.], [0., None, 0.], 'missing_data'),
                                 ([0., 1., 2.], [0., float('nan'), 0.], 'nonfinite'),
                                 ([0., 1., 2.], [0., float('inf'), 0.], 'nonfinite'),
                                 ([0., 0., 2.], [0., 2., 0.], 'non_increasing_time'),
                                 ([0., 1., 10.], [0., 2., 0.], 'acquisition_gap')]:
            with self.subTest(code=code), self.assertRaisesRegex(ValueError, code):
                self.api.integrate_peak(rt, signal, bounds_min=(0., rt[-1]), baseline='none', max_gap_min=2.)

    def test_peak_selection_is_rt_constrained_and_ambiguous_peaks_require_review(self):
        peaks = [dict(apex_rt_min=1., area=1000.), dict(apex_rt_min=5., area=10.)]
        result = self.api.select_peak(peaks, expected_rt_min=5., tolerance_min=0.2)
        self.assertEqual(result['peak']['area'], 10.)
        result = self.api.select_peak(peaks + [dict(apex_rt_min=5.1, area=12.)], expected_rt_min=5., tolerance_min=0.2)
        self.assertEqual(result['status'], 'ambiguous_peak')
        self.assertIsNone(result['peak'])

    def test_fit_exact_line_with_free_intercept(self):
        result = self.api.fit_calibration(standards(), weighting='none', intercept='free')
        self.assertAlmostEqual(result['slope'], 2.)
        self.assertAlmostEqual(result['intercept'], 1.)
        self.assertEqual(result['range'], [1., 4.])
        self.assertEqual(result['n_levels'], 3)
        self.assertEqual(result['n_points'], 3)
        self.assertTrue(all(abs(x) < 1e-10 for x in result['backcalc_bias_pct']))

    def test_weighted_fits_have_independent_numeric_oracles(self):
        rows = standards()
        for row, response in zip(rows, [3., 5., 10.]):
            row['response'] = response
        for weighting, slope, intercept in [('none', 33/14, 1/2), ('1/x', 30/13, 8/13), ('1/x2', 9/4, 5/7)]:
            with self.subTest(weighting=weighting):
                result = self.api.fit_calibration(rows, weighting=weighting, intercept='free')
                self.assertAlmostEqual(result['slope'], slope)
                self.assertAlmostEqual(result['intercept'], intercept)

    def test_forced_origin_requires_explicit_option(self):
        result = self.api.fit_calibration(standards(), weighting='none', intercept='zero')
        self.assertEqual(result['intercept'], 0.)
        self.assertAlmostEqual(result['slope'], 7/3)

    def test_bad_calibration_inputs_are_rejected(self):
        cases = [('concentration', None, 'missing_concentration'),
                 ('concentration', float('nan'), 'nonfinite'),
                 ('concentration', -1., 'negative_concentration'),
                 ('response', float('inf'), 'nonfinite'),
                 ('batch_id', 'B2', 'batch_mismatch'),
                 ('method_fingerprint', 'M2', 'method_mismatch'),
                 ('concentration_unit', 'mg/L', 'unit_mismatch'),
                 ('analyte_id', 'B', 'analyte_mismatch')]
        for key, value, code in cases:
            rows = standards(); rows[0][key] = value
            with self.subTest(code=code), self.assertRaisesRegex(ValueError, code):
                self.api.fit_calibration(rows, weighting='none', intercept='free')

    def test_zero_calibration_level_is_rejected_for_inverse_weights(self):
        for weighting in ['1/x', '1/x2']:
            rows = standards(); rows[0]['concentration'] = 0.
            with self.subTest(weighting=weighting), self.assertRaisesRegex(ValueError, 'zero_weight_domain'):
                self.api.fit_calibration(rows, weighting=weighting, intercept='free')

    def test_blank_and_qc_are_not_fit_points(self):
        rows = standards()
        for role in ['blank', 'qc', 'unknown']:
            extra = copy.deepcopy(rows[0]); extra.update(role=role, run_id=role, response=1e9)
            rows.append(extra)
        result = self.api.fit_calibration(rows, weighting='none', intercept='free')
        self.assertAlmostEqual(result['slope'], 2.)
        self.assertEqual(result['n_points'], 3)

    def test_replicates_do_not_inflate_number_of_levels(self):
        rows = standards()
        extra = copy.deepcopy(rows[0]); extra['run_id'] = 'replicate'
        result = self.api.fit_calibration(rows + [extra], weighting='none', intercept='free')
        self.assertEqual(result['n_points'], 4)
        self.assertEqual(result['n_levels'], 3)
        for row in rows:
            row['concentration'] = 1.
        with self.assertRaisesRegex(ValueError, 'insufficient_levels'):
            self.api.fit_calibration(rows, weighting='none', intercept='free')

    def test_outliers_are_not_silently_deleted(self):
        rows = standards(); rows[-1]['response'] = 1000.
        result = self.api.fit_calibration(rows, weighting='none', intercept='free')
        self.assertEqual(result['n_points'], 3)
        self.assertEqual(result['excluded_run_ids'], [])
        rows[-1]['include'] = False
        with self.assertRaisesRegex(ValueError, 'exclusion_reason_required'):
            self.api.fit_calibration(rows, weighting='none', intercept='free')

    def test_internal_standard_uses_area_ratio_and_requires_positive_area(self):
        self.assertAlmostEqual(self.api.response_value(100., mode='internal', is_area=20.), 5.)
        self.assertEqual(self.api.response_value(100., mode='external'), 100.)
        for value in [None, 0., -1.]:
            with self.subTest(is_area=value), self.assertRaisesRegex(ValueError, 'invalid_internal_standard'):
                self.api.response_value(100., mode='internal', is_area=value)

    def test_backcalculation_applies_dilution_once(self):
        model = self.api.fit_calibration(standards(), weighting='none', intercept='free')
        result = self.api.quantify(5., model, dilution_factor=10., batch_id='B1', method_fingerprint='M1')
        self.assertAlmostEqual(result['vial_concentration'], 2.)
        self.assertAlmostEqual(result['original_concentration'], 20.)
        self.assertEqual(result['status'], 'ok')

    def test_out_of_range_is_flagged_and_not_reported_as_valid(self):
        model = self.api.fit_calibration(standards(), weighting='none', intercept='free')
        for response, code in [(11., 'above_calibration_range'), (2., 'below_calibration_range'), (0., 'negative_backcalc')]:
            with self.subTest(code=code):
                result = self.api.quantify(response, model, dilution_factor=1., batch_id='B1', method_fingerprint='M1')
                self.assertEqual(result['status'], code)
                self.assertIsNone(result['original_concentration'])

    def test_quantification_rejects_other_batch_method_and_invalid_dilution(self):
        model = self.api.fit_calibration(standards(), weighting='none', intercept='free')
        for batch, method, dilution, code in [('B2', 'M1', 1., 'batch_mismatch'),
                                               ('B1', 'M2', 1., 'method_mismatch'),
                                               ('B1', 'M1', 0., 'invalid_dilution')]:
            with self.subTest(code=code), self.assertRaisesRegex(ValueError, code):
                self.api.quantify(5., model, dilution_factor=dilution, batch_id=batch, method_fingerprint=method)

    def test_qualifier_ratio_is_relative_fraction_and_requires_coelution(self):
        for qual, delta_rt, status in [(40., 0.01, 'ok'), (60., 0.01, 'ion_ratio_fail'), (40., 0.2, 'rt_mismatch')]:
            result = self.api.check_qualifier(100., qual, reference_ratio=0.4,
                                              relative_tolerance=0.2, apex_delta_min=delta_rt, rt_tolerance_min=0.05)
            self.assertEqual(result['status'], status)
        result = self.api.check_qualifier(0., 0., reference_ratio=0.4,
                                          relative_tolerance=0.2, apex_delta_min=0., rt_tolerance_min=0.05)
        self.assertEqual(result['status'], 'invalid_quantifier')

    def test_output_may_not_be_inside_raw_or_existing_result_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); raw = root / 'run.d'; raw.mkdir()
            existing = root / 'previous_results'; existing.mkdir()
            for output in [raw, raw / 'results', existing]:
                with self.subTest(output=output), self.assertRaisesRegex(ValueError, 'unsafe_output'):
                    self.api.validate_output_path(output, raw_datasets=[raw])
            self.api.validate_output_path(root / 'new_results', raw_datasets=[raw])


if __name__ == '__main__':
    unittest.main()
