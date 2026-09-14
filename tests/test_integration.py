"""Integration behaviour beyond the shared contracts, plus the real standards."""
import unittest

import context

from mrm_quant import integration, reader_adapter, transitions

# Configuration used by these tests for the September batch. Chosen from the
# inspected traces of this batch, not a general recommendation.
EXPECTED_RT_MIN = 14.43
RT_TOLERANCE_MIN = 0.08
SEARCH_MARGIN_MIN = 0.25
FIND_OPTIONS = dict(min_relative_height=0.05, min_separation_min=0.1)


class IntegrateWindowTests(unittest.TestCase):
    def test_flat_signal_with_endpoint_baseline_is_not_a_quantification_peak(self):
        result = integration.integrate_peak([0., 1., 2.], [5., 5., 5.],
                                            bounds_min=(0., 2.),
                                            baseline='linear_endpoints', max_gap_min=2.)
        self.assertEqual(result['status'], 'non_positive_area')
        self.assertAlmostEqual(result['area'], 0.)

    def test_window_outside_the_acquired_data_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'window_outside_data'):
            integration.integrate_peak([1., 2., 3.], [1., 2., 1.],
                                       bounds_min=(0.5, 3.), baseline='none', max_gap_min=2.)

    def test_unknown_baseline_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'invalid_baseline'):
            integration.integrate_peak([0., 1., 2.], [0., 2., 0.],
                                       bounds_min=(0., 2.), baseline='rolling', max_gap_min=2.)

    def test_too_few_points_do_not_produce_an_area(self):
        result = integration.integrate_peak([0., 1., 2.], [0., 2., 0.],
                                            bounds_min=(0., 0.5), baseline='none',
                                            max_gap_min=2.)
        self.assertEqual(result['status'], 'insufficient_points')
        self.assertIsNone(result['area'])

    def test_no_candidate_in_the_window_is_reported(self):
        peaks = [dict(apex_rt_min=1., area=10.)]
        result = integration.select_peak(peaks, expected_rt_min=5., tolerance_min=0.2)
        self.assertEqual(result['status'], 'no_peak')
        self.assertIsNone(result['peak'])

    def test_area_is_signed_after_baseline_subtraction(self):
        # A dip below the endpoint baseline must not be clipped to zero.
        result = integration.integrate_peak([0., 1., 2.], [10., 6., 10.],
                                            bounds_min=(0., 2.),
                                            baseline='linear_endpoints', max_gap_min=2.)
        self.assertAlmostEqual(result['area'], -4.)
        self.assertEqual(result['status'], 'non_positive_area')


class SeptemberStandardsTests(unittest.TestCase):
    def trace(self, name):
        path = context.dataset_path('TSJ-0907/LQ/%s.d' % name)
        if not path.is_dir():
            self.skipTest('Local raw dataset unavailable: ' + name)
        with reader_adapter.RunReader(path) as run:
            channel = next(c for c in run.channels
                           if c['channel_type'] == 'mrm' and c['product_mz'] == 93.)
        return transitions.extract_run(path, transitions.selector_from_channel(channel))

    def area(self, name):
        trace = self.trace(name)
        peaks = integration.find_peaks(
            trace['rt_min'], trace['intensity'],
            rt_range=(EXPECTED_RT_MIN - SEARCH_MARGIN_MIN,
                      EXPECTED_RT_MIN + SEARCH_MARGIN_MIN), **FIND_OPTIONS)
        chosen = integration.select_peak(peaks, expected_rt_min=EXPECTED_RT_MIN,
                                         tolerance_min=RT_TOLERANCE_MIN)
        return chosen

    def test_every_standard_yields_exactly_one_peak_in_the_window(self):
        for level in range(1, 11):
            with self.subTest(level=level):
                chosen = self.area('STD_S%d' % level)
                self.assertEqual(chosen['status'], 'ok')
                self.assertEqual(chosen['candidates'], 1)
                self.assertEqual(chosen['peak']['status'], 'ok')
                self.assertGreater(chosen['peak']['area'], 0.)
                self.assertEqual(chosen['peak']['area_unit'], 'stored_intensity*min')

    def test_areas_decrease_along_the_two_fold_dilution_series(self):
        areas = [self.area('STD_S%d' % level)['peak']['area'] for level in range(1, 11)]
        self.assertEqual(areas, sorted(areas, reverse=True))

    def test_samples_with_the_same_peak_are_integrated_the_same_way(self):
        for name in ('S21', 'S22', 'S23', 'S24'):
            with self.subTest(run=name):
                chosen = self.area(name)
                self.assertEqual(chosen['status'], 'ok')
                self.assertGreater(chosen['peak']['area'], 0.)


if __name__ == '__main__':
    unittest.main()
