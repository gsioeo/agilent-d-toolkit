"""Transition extraction beyond the shared contracts: grouping and real runs."""
import unittest

import context

from mrm_quant import reader_adapter, transitions


def mrm_record(scan_id, rt, *, method_id=1, segment=1, intensities=(1., 2., 10.)):
    return dict(scan_id=scan_id, rt_min=rt, scan_type=256, ms_level=2,
                precursor_mz=204., collision_energy_ev=30., polarity=0,
                time_segment_id=segment, scan_method_id=method_id,
                cycle_number=scan_id, product_mz=[68., 81., 93.],
                intensity=list(intensities), declared_product_mz=[68., 81., 93.])


class GroupingTests(unittest.TestCase):
    def test_two_mrm_groups_are_not_summed_implicitly(self):
        records = [mrm_record(1, 1.), mrm_record(2, 1.1, method_id=3)]
        with self.assertRaisesRegex(ValueError, 'multiple_acquisition_groups'):
            transitions.separate_tic(records)

    def test_selector_from_channel_round_trips(self):
        channel = dict(channel_type='mrm', channel_id='c', precursor_mz=204.,
                       product_mz=93., collision_energy_ev=30., polarity=0,
                       time_segment_id=1, scan_method_id=1)
        selector = transitions.selector_from_channel(channel, product_tol_da=0.02)
        result = transitions.extract_transition([mrm_record(1, 1.)], selector)
        self.assertEqual(result['intensity'], [10.])
        self.assertEqual(result['selector']['product_tol_da'], 0.02)

    def test_ms1_channel_cannot_be_used_as_a_transition(self):
        channel = dict(channel_type='ms1', channel_id='ms1', time_segment_id=1,
                       scan_method_id=2, mz_low=50., mz_high=500.)
        with self.assertRaisesRegex(ValueError, 'invalid_selector'):
            transitions.selector_from_channel(channel)


class RealRunTests(unittest.TestCase):
    def setUp(self):
        self.path = context.dataset_path('TSJ-0907/LQ/STD_S1.d')
        if not self.path.is_dir():
            self.skipTest('Local raw dataset unavailable')

    def selector(self, product):
        with reader_adapter.RunReader(self.path) as run:
            channel = next(c for c in run.channels
                           if c['channel_type'] == 'mrm' and c['product_mz'] == product)
            return transitions.selector_from_channel(channel)

    def test_every_declared_channel_extracts_native_sampling_only(self):
        for product in (93., 81., 68.):
            with self.subTest(product=product):
                result = transitions.extract_run(self.path, self.selector(product))
                self.assertEqual(result['point_count'], 2454)
                self.assertEqual(set(result['point_status']), {'ok'})
                self.assertTrue(all(rt >= 10. for rt in result['rt_min']))

    def test_unacquired_product_is_refused_for_a_real_run(self):
        selector = self.selector(93.)
        selector['product_mz'] = 189.
        with self.assertRaisesRegex(ValueError, 'transition_not_acquired'):
            transitions.extract_run(self.path, selector)

    def test_mrm_sum_equals_the_three_transitions_and_ms1_stays_separate(self):
        separated = transitions.separate_run(self.path)
        self.assertEqual(separated['mrm_sum']['point_count'], 2454)
        self.assertEqual(separated['ms1_tic']['point_count'], 2454)
        self.assertEqual(separated['mrm_sum']['group'], (1, 1))
        self.assertEqual(separated['ms1_tic']['group'], (1, 2))
        traces = [transitions.extract_run(self.path, self.selector(p))['intensity']
                  for p in (93., 81., 68.)]
        for index in (0, 1000, 2453):
            self.assertAlmostEqual(separated['mrm_sum']['intensity'][index],
                                   sum(trace[index] for trace in traces), places=6)


if __name__ == '__main__':
    unittest.main()
