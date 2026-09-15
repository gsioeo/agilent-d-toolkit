"""Calibration behaviour beyond the shared synthetic contracts."""
import unittest

from mrm_quant import calibration


def rows(responses, concentrations=(1., 2., 4.)):
    return [{'run_id': 'cal%d' % index, 'batch_id': 'B1', 'analyte_id': 'A',
             'method_fingerprint': 'M1', 'role': 'calibration',
             'concentration': concentration, 'concentration_unit': 'ng/mL',
             'response': response, 'include': True, 'exclusion_reason': None}
            for index, (concentration, response)
            in enumerate(zip(concentrations, responses))]


class FitDiagnosticsTests(unittest.TestCase):
    def test_equation_and_r2_match_an_independent_calculation(self):
        model = calibration.fit_calibration(rows([3., 5., 10.]), weighting='none',
                                            intercept='free')
        self.assertEqual(model['equation'], 'A = 2.35714 * C + 0.5')
        # Least squares on x=[1,2,4], y=[3,5,10]: residuals 1/7, -3/14, 1/14 and
        # total sum of squares 26, so r2 = 1 - (1/2)/(14*26) * 14 = 0.9972527...
        self.assertAlmostEqual(model['r2'], 1 - (1/7)**2 * (1 + 2.25 + 0.25) / 26, places=12)
        self.assertAlmostEqual(model['r2'], 0.99725274725274725, places=12)
        self.assertAlmostEqual(model['r2_weighted'], model['r2'], places=12)
        self.assertEqual(model['weight_definition'], 'w = 1')
        self.assertEqual(model['model'], 'linear')

    def test_perfect_line_has_unit_r2(self):
        model = calibration.fit_calibration(rows([3., 5., 9.]), weighting='none',
                                            intercept='free')
        self.assertAlmostEqual(model['r2'], 1.)
        self.assertAlmostEqual(model['slope'], 2.)

    def test_weighted_r2_uses_the_fit_weights(self):
        table = rows([3., 5., 10.])
        model = calibration.fit_calibration(table, weighting='1/x', intercept='free')
        weights = [point['weight'] for point in model['points']]
        responses = [point['response'] for point in model['points']]
        predicted = [model['slope'] * point['concentration'] + model['intercept']
                     for point in model['points']]
        mean = sum(w * y for w, y in zip(weights, responses)) / sum(weights)
        expected = 1 - (sum(w * (y - p) ** 2 for w, y, p in zip(weights, responses, predicted))
                        / sum(w * (y - mean) ** 2 for w, y in zip(weights, responses)))
        self.assertAlmostEqual(model['r2_weighted'], expected, places=12)
        self.assertNotAlmostEqual(model['r2_weighted'], model['r2'])

    def test_falling_response_is_not_a_calibration(self):
        with self.assertRaisesRegex(ValueError, 'non_positive_slope'):
            calibration.fit_calibration(rows([10., 5., 1.]), weighting='none',
                                        intercept='free')

    def test_unknown_weighting_and_intercept_are_refused(self):
        with self.assertRaisesRegex(ValueError, 'invalid_weighting'):
            calibration.fit_calibration(rows([3., 5., 10.]), weighting='1/y',
                                        intercept='free')
        with self.assertRaisesRegex(ValueError, 'invalid_intercept'):
            calibration.fit_calibration(rows([3., 5., 10.]), weighting='none',
                                        intercept='through_origin')

    def test_excluded_point_is_recorded_with_its_reason(self):
        table = rows([3., 5., 10., 11.], concentrations=(1., 2., 4., 8.))
        table[-1].update(include=False, exclusion_reason='syringe misfire')
        model = calibration.fit_calibration(table, weighting='none', intercept='free')
        self.assertEqual(model['excluded_run_ids'], ['cal3'])
        self.assertEqual(model['exclusions'][0]['reason'], 'syringe misfire')
        self.assertEqual(model['n_points'], 3)

    def test_internal_standard_response_needs_a_positive_area(self):
        self.assertAlmostEqual(calibration.response_value(50., mode='internal',
                                                          is_area=25.), 2.)
        with self.assertRaisesRegex(ValueError, 'invalid_response_mode'):
            calibration.response_value(50., mode='ratio')


if __name__ == '__main__':
    unittest.main()
