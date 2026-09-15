"""v2 MRM search, MassQL mapping and independent-reader comparison."""
import unittest

import context

from mrm_quant import massql_adapter, reader_validation, search


def channel(**changes):
    value = {'channel_type': 'mrm', 'channel_id': 'seg1_m1_204to93_ce30',
             'compound_name': 'target', 'precursor_mz': 204., 'product_mz': 93.,
             'collision_energy_ev': 30., 'polarity': 0, 'time_segment_id': 1,
             'scan_method_id': 1, 'frame_count': 20}
    value.update(changes)
    return value


def record(scan_id, ms_level, mz, intensity, **changes):
    value = {'scan_id': scan_id, 'rt_min': float(scan_id), 'scan_type': 256 if ms_level == 2 else 1,
             'ms_level': ms_level, 'precursor_mz': 204. if ms_level == 2 else 0.,
             'collision_energy_ev': 30. if ms_level == 2 else 0., 'polarity': 0,
             'time_segment_id': 1, 'scan_method_id': ms_level, 'cycle_number': scan_id,
             'product_mz': list(mz), 'intensity': list(intensity),
             'declared_product_mz': list(mz), 'point_count': len(mz)}
    value.update(changes)
    return value


class ChannelSearchTests(unittest.TestCase):
    def test_optional_criteria_compose_without_rounding(self):
        channels = [channel(), channel(channel_id='other', product_mz=81.),
                    channel(channel_id='ms1', channel_type='ms1')]
        criteria = {'precursor_mz': 204.004, 'product_mz': 93.004,
                    'collision_energy_ev': 30.005, 'polarity': 0,
                    'time_segment_id': 1, 'scan_method_id': 1,
                    'precursor_tol_da': .005, 'product_tol_da': .005,
                    'ce_tol_ev': .01}
        self.assertEqual([item['channel_id'] for item in
                          search.search_channels(channels, criteria)],
                         ['seg1_m1_204to93_ce30'])

    def test_multiple_matches_remain_multiple(self):
        channels = [channel(), channel(channel_id='second', product_mz=81.)]
        self.assertEqual(len(search.search_channels(channels, {'precursor_mz': 204.,
                                                                'precursor_tol_da': .01})), 2)


class FakeFrame:
    def __init__(self, rows, columns):
        self.rows = rows
        self.columns = columns


class FakePandas:
    DataFrame = FakeFrame


class MassQLMappingTests(unittest.TestCase):
    def test_interleaved_scans_map_to_separate_frames(self):
        records = [record(1, 1, [50., 51.], [2., 6.]),
                   record(2, 2, [81., 93.], [0., 10.])]
        ms1, ms2 = massql_adapter.records_to_dataframes(records, FakePandas)
        self.assertEqual(len(ms1.rows), 2)
        self.assertEqual(len(ms2.rows), 2)
        self.assertEqual(ms2.rows[1]['precmz'], 204.)
        self.assertEqual(ms2.rows[1]['ms1scan'], 1)
        self.assertEqual(ms2.rows[1]['i_norm'], 1.)
        self.assertEqual(ms2.rows[0]['i_tic_norm'], 0.)

    def test_vendor_context_is_not_lost(self):
        rows = massql_adapter.scan_context('run', [record(2, 2, [93.], [5.])])
        self.assertEqual(rows[0]['time_segment_id'], 1)
        self.assertEqual(rows[0]['scan_method_id'], 2)
        self.assertEqual(rows[0]['collision_energy_ev'], 30.)


class FakeCentroid:
    def __init__(self, times, scans):
        self.xlabels = times
        self.scans = scans

    def scan(self, index):
        return self.scans[index]


class ReaderComparisonTests(unittest.TestCase):
    def test_equal_scans_pass_within_tolerance(self):
        records = [record(1, 1, [50., 51.], [2., 6.])]
        independent = FakeCentroid([1.0000001], [([50., 51.], [2., 6.000001])])
        rows = reader_validation.compare_records(
            records, independent, rt_tolerance=1e-6, intensity_rtol=1e-6)
        self.assertEqual(rows[0]['status'], 'ok')

    def test_scan_and_peak_mismatches_are_explicit(self):
        records = [record(1, 1, [50., 51.], [2., 6.]),
                   record(2, 2, [93.], [5.])]
        independent = FakeCentroid([1.], [([50.], [2.])])
        rows = reader_validation.compare_records(records, independent)
        self.assertEqual([row['status'] for row in rows],
                         ['peak_count_mismatch', 'missing_rainbow'])

    def test_integer_only_independent_values_are_classified(self):
        records = [record(1, 1, [50., 51.], [2.25, 6.75])]
        independent = FakeCentroid([1.], [([50., 51.], [2, 6])])
        rows = reader_validation.compare_records(records, independent)
        self.assertEqual(rows[0]['status'], 'integer_quantized')


class OptionalIntegrationTests(unittest.TestCase):
    def test_massql_engine_accepts_direct_dataframes(self):
        try:
            from massql import msql_engine
        except ImportError:
            self.skipTest('MassQL unavailable')
        records = [record(1, 1, [50., 93.], [1., 10.]),
                   record(2, 2, [81., 93.], [2., 20.])]
        ms1, ms2 = massql_adapter.records_to_dataframes(records)
        result = msql_engine.process_query(
            'QUERY scaninfo(MS2DATA) WHERE MS2PREC=204:TOLERANCEMZ=0.01 '
            'AND MS2PROD=93:TOLERANCEMZ=0.01',
            'synthetic.mzML', ms1_df=ms1, ms2_df=ms2)
        self.assertEqual(result['scan'].tolist(), [2])
        rejected = msql_engine.process_query(
            'QUERY scaninfo(MS2DATA) WHERE MS2PROD=81:TOLERANCEMZ=0.01:'
            'INTENSITYPERCENT=20',
            'synthetic.mzML', ms1_df=ms1, ms2_df=ms2)
        self.assertEqual(len(rejected), 0)

    def test_real_run_agrees_with_rainbow(self):
        try:
            import rainbow
        except ImportError:
            self.skipTest('rainbow unavailable')
        dataset = context.dataset_path('mrm_primary')
        if not dataset.is_dir():
            self.skipTest('Local raw dataset unavailable')
        from mrm_quant import reader_adapter
        with reader_adapter.RunReader(dataset) as run:
            records = run.records()
        directory = rainbow.read(str(dataset), centroid=True, display_precision=8)
        centroid = next(item for item in directory.datafiles
                        if item.name.lower() == 'mspeak.bin')
        rows = reader_validation.compare_records(records, centroid)
        self.assertEqual(len(rows), 4908)
        self.assertEqual({row['status'] for row in rows},
                         {'ok', 'integer_quantized'})


if __name__ == '__main__':
    unittest.main()
