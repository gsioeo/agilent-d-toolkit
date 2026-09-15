"""Read-only characterization tests; no new extraction implementation."""
import hashlib
import json
from pathlib import Path
import unittest

import context

REPO = context.REPO
legacy = context.load_legacy('agilent_d')


class ExistingReaderTests(unittest.TestCase):
    def dataset(self, key):
        path = context.dataset_path(key)
        if not path.is_dir():
            self.skipTest('Local raw fixture unavailable: ' + key)
        ds = legacy.AgilentDotD(str(path))
        self.addCleanup(ds.close)
        return ds

    def test_original_scripts_unchanged(self):
        hashes = json.loads((REPO / 'tests/legacy_sha256.json').read_text())
        for name, expected in hashes.items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((REPO / name).read_bytes()).hexdigest(), expected)

    def test_full_scan_fixture_contains_only_ms1_records(self):
        ds = self.dataset('ms1_reference')
        self.assertEqual({(s.scan_type, s.ms_level) for s in ds.scans}, {(1, 1)})

    def test_mrm_fixture_interleaves_two_modes_in_each_cycle(self):
        ds = self.dataset('mrm_primary')
        self.assertGreater(len(ds.scans), 0)
        self.assertEqual(len(ds.scans) % 2, 0)
        for mrm, ms1 in zip(ds.scans[::2], ds.scans[1::2]):
            self.assertEqual((mrm.scan_type, mrm.ms_level), (256, 2))
            self.assertEqual((ms1.scan_type, ms1.ms_level), (1, 1))
            self.assertEqual(mrm.cycle_number, ms1.cycle_number)
            self.assertLess(mrm.scan_time, ms1.scan_time)

    def test_binary_blocks_are_structurally_consistent(self):
        ds = self.dataset('mrm_primary')
        size = (Path(ds.acqdata) / 'MSPeak.bin').stat().st_size
        for scan, mz, abundance in ds.iter_spectra():
            self.assertEqual(scan.byte_count, scan.point_count * 8)
            self.assertLessEqual(scan.spectrum_offset + scan.byte_count, size)
            self.assertEqual(tuple(sorted(mz)), tuple(mz))
            self.assertAlmostEqual(sum(abundance), scan.tic,
                                   delta=max(1e-6, abs(scan.tic) * 1e-10))

    def test_eic_and_xic_paths_agree(self):
        try:
            import numpy  # noqa: F401
        except ImportError:
            self.skipTest('numpy required by existing eic()')
        ds = self.dataset('mrm_primary')
        mz, abundance = ds.spectrum(0)
        target = mz[-1]
        rt, signals = ds.eic([target], tol=0.3)
        lo, hi = ds.eic_window(target, 0.3)
        slow_rt, slow = ds.xic(lo, hi)
        self.assertEqual(list(rt), slow_rt)
        self.assertEqual(list(signals[0]), slow)


if __name__ == '__main__':
    unittest.main()
