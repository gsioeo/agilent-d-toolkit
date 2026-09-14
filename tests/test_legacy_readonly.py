"""Read-only characterization tests; no new extraction implementation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

import context

REPO = context.REPO
legacy = context.load_legacy('agilent_d')


class ExistingReaderTests(unittest.TestCase):
    def dataset(self, relative):
        path = context.dataset_path(relative)
        if not path.is_dir():
            self.skipTest('Local raw dataset unavailable: ' + relative)
        ds = legacy.AgilentDotD(str(path))
        self.addCleanup(ds.close)
        return ds

    def test_original_scripts_unchanged(self):
        hashes = json.loads((REPO / 'tests/legacy_sha256.json').read_text())
        for name, expected in hashes.items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((REPO / name).read_bytes()).hexdigest(), expected)

    def test_july_is_ms1_even_though_instrument_is_qqq(self):
        ds = self.dataset('20260729/STD.d')
        self.assertEqual(ds.ms_method['msInstrument'], 'QQQ')
        self.assertEqual({(s.scan_type, s.ms_level) for s in ds.scans}, {(1, 1)})
        self.assertEqual(len(ds.scans), 6269)

    def test_september_interleaves_two_modes_in_each_cycle(self):
        ds = self.dataset('TSJ-0907/LQ/STD_S1.d')
        self.assertEqual(len(ds.scans), 4908)
        for mrm, ms1 in zip(ds.scans[::2], ds.scans[1::2]):
            self.assertEqual((mrm.scan_type, mrm.ms_level, mrm.scan_method_id), (256, 2, 1))
            self.assertEqual((ms1.scan_type, ms1.ms_level, ms1.scan_method_id), (1, 1, 2))
            self.assertEqual(mrm.cycle_number, ms1.cycle_number)
            self.assertLess(mrm.scan_time, ms1.scan_time)

    def test_method_transitions_differ_between_batches(self):
        for relative, products in [('20260819-tsj/TSJ-0819-1.d', {69., 93., 189.}),
                                   ('TSJ-0907/LQ/STD_S1.d', {68., 81., 93.})]:
            with self.subTest(batch=relative):
                ds = self.dataset(relative)
                elements = ds.ms_method['timeSegments'][0]['scanSegments'][0]['scanElements']
                self.assertEqual({float(e['ms1LowMz']) for e in elements}, {204.})
                self.assertEqual({float(e['ms2LowMz']) for e in elements}, products)
                self.assertEqual({float(e['collisionEnergy']) for e in elements}, {30.})
                self.assertEqual({float(e['dwell']) for e in elements}, {100.})
                self.assertEqual({e['isISTD'] for e in elements}, {'false'})

    def test_binary_product_order_is_not_method_element_order(self):
        ds = self.dataset('TSJ-0907/LQ/STD_S1.d')
        elements = ds.ms_method['timeSegments'][0]['scanSegments'][0]['scanElements']
        self.assertEqual([float(e['ms2LowMz']) for e in elements], [93., 81., 68.])
        mz, abundance = ds.spectrum(0)
        self.assertEqual(mz, (68., 81., 93.))
        self.assertEqual(abundance, (0.11013031005859375, 5.686073303222656, 1121.11181640625))
        self.assertEqual(ds.scans[0].mz_of_interest, 204.)

    def test_representative_binary_blocks_and_tic_consistency(self):
        # Structural consistency supports this local layout; it does not prove
        # intensity units or equivalence to a vendor chromatogram export.
        for relative in ['20260729/STD.d',
                         '20260819-tsj/TSJ-0819-1.d', 'TSJ-0907/LQ/STD_S1.d']:
            with self.subTest(dataset=relative):
                ds = self.dataset(relative)
                size = (Path(ds.acqdata) / 'MSPeak.bin').stat().st_size
                for s, mz, abundance in ds.iter_spectra():
                    self.assertEqual(s.byte_count, s.point_count * 8)
                    self.assertGreaterEqual(s.spectrum_offset, 0)
                    self.assertLessEqual(s.spectrum_offset + s.byte_count, size)
                    self.assertEqual(tuple(sorted(mz)), tuple(mz))
                    self.assertAlmostEqual(sum(abundance), s.tic, delta=max(1e-6, abs(s.tic) * 1e-10))
                    if s.scan_type == 256:
                        self.assertEqual(s.spectrum_format_id, 3)
                        self.assertEqual(s.point_count, 3)
                        self.assertEqual(s.mz_of_interest, 204.)

    def test_tic_returns_every_record_without_mode_separation(self):
        ds = self.dataset('TSJ-0907/LQ/STD_S1.d')
        rt, tic = ds.tic()
        self.assertEqual(rt, [s.scan_time for s in ds.scans])
        self.assertEqual(tic, [s.tic for s in ds.scans])
        self.assertEqual(len(tic), 4908)  # Only 2454 of these are MS1 TIC points.

    def test_existing_eic_contains_ms1_signal_and_is_not_a_transition(self):
        if importlib.util.find_spec('numpy') is None:
            self.skipTest('numpy required by existing eic()')
        ds = self.dataset('TSJ-0907/LQ/STD_S1.d')
        rt, signals = ds.eic([93.], tol=0.3)
        self.assertEqual(len(rt), 4908)
        self.assertEqual(signals.shape, (1, 4908))
        self.assertTrue(any(signals[0, i] > 0 for i, s in enumerate(ds.scans) if s.ms_level == 1))
        self.assertEqual(signals[0, 0], 1121.11181640625)
        lo, hi = ds.eic_window(93., 0.3)
        slow_rt, slow = ds.xic(lo, hi)
        self.assertEqual(list(rt), slow_rt)
        self.assertEqual(list(signals[0]), slow)

    def test_all_september_standards_have_no_declared_level(self):
        for i in range(1, 11):
            with self.subTest(standard=i):
                ds = self.dataset(f'TSJ-0907/LQ/STD_S{i}.d')
                self.assertEqual(ds.sample_info.get('Level Name'), '')
                sequence = ET.parse(Path(ds.acqdata) / 'sequence.xml').getroot().find('Sequence')
                self.assertEqual(sequence.findtext('LevelName'), '')
                self.assertEqual(sequence.findtext('SampleType'), 'Sample')


if __name__ == '__main__':
    unittest.main()
