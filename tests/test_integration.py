"""Integration behaviour beyond the shared synthetic contracts."""
import unittest

from mrm_quant import integration


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


if __name__ == '__main__':
    unittest.main()
