"""Safety checks for regenerated legacy-ingestion artifacts."""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from ingest.output import commit_output, discard_output, begin_output, unique_dataset_basenames
from ingest import eic


class LegacyOutputTests(unittest.TestCase):
    def test_duplicate_dot_d_basenames_are_rejected_case_insensitively(self):
        with self.assertRaisesRegex(ValueError, 'duplicate_dataset_basename'):
            unique_dataset_basenames(['/a/sample.d', '/b/SAMPLE.D'])

    def test_commit_replaces_generated_files_and_preserves_user_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / 'ingested'
            out.mkdir()
            user_file = out / 'notes.txt'
            user_file.write_text('keep', encoding='utf-8')
            stage = begin_output(out)
            (stage / 'manifest.json').write_text('new', encoding='utf-8')
            commit_output(stage, out, group='ingest')
            self.assertEqual((out / 'manifest.json').read_text(encoding='utf-8'), 'new')
            self.assertEqual(user_file.read_text(encoding='utf-8'), 'keep')

    def test_discarded_stage_does_not_touch_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / 'ingested'
            out.mkdir()
            stage = begin_output(out)
            (stage / 'partial.csv').write_text('partial', encoding='utf-8')
            discard_output(stage)
            self.assertFalse((out / 'partial.csv').exists())

    def test_known_stale_artifacts_are_removed_but_nested_user_files_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / 'ingested'
            (out / 'data' / 'old').mkdir(parents=True)
            (out / 'mzml').mkdir()
            (out / 'notes').mkdir()
            (out / 'data' / 'old' / 'tic.csv').write_text('old', encoding='utf-8')
            (out / 'data' / 'old' / 'user.txt').write_text('keep', encoding='utf-8')
            (out / 'mzml' / 'old.mzML').write_text('old', encoding='utf-8')
            (out / 'notes' / 'user.txt').write_text('keep', encoding='utf-8')
            stage = begin_output(out)
            (stage / 'data' / 'new').mkdir(parents=True)
            (stage / 'data' / 'new' / 'tic.csv').write_text('new', encoding='utf-8')
            commit_output(stage, out, group='ingest')
            self.assertFalse((out / 'data' / 'old' / 'tic.csv').exists())
            self.assertFalse((out / 'mzml' / 'old.mzML').exists())
            self.assertEqual((out / 'data' / 'old' / 'user.txt').read_text(encoding='utf-8'), 'keep')
            self.assertEqual((out / 'notes' / 'user.txt').read_text(encoding='utf-8'), 'keep')

    def test_tracked_rerun_preserves_untracked_file_with_generated_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'ingested'
            stage = begin_output(out)
            (stage / 'data' / 'first').mkdir(parents=True)
            (stage / 'data' / 'first' / 'tic.csv').write_text('first', encoding='utf-8')
            commit_output(stage, out, group='ingest')
            user_file = out / 'data' / 'user' / 'tic.csv'
            user_file.parent.mkdir()
            user_file.write_text('keep', encoding='utf-8')
            stage = begin_output(out)
            (stage / 'data' / 'second').mkdir(parents=True)
            (stage / 'data' / 'second' / 'tic.csv').write_text('second', encoding='utf-8')
            commit_output(stage, out, group='ingest')
            self.assertFalse((out / 'data' / 'first' / 'tic.csv').exists())
            self.assertEqual(user_file.read_text(encoding='utf-8'), 'keep')

    def test_eic_discards_stage_on_unexpected_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / 'eic'
            with mock.patch.object(eic, 'find_datasets', return_value=['/data/run.d']), \
                    mock.patch.object(eic, 'AgilentDotD', side_effect=RuntimeError('broken')):
                with self.assertRaisesRegex(RuntimeError, 'broken'):
                    eic.main(['--source', str(root), '--out', str(out),
                              '--mz', '1', '--no-plots'])
            self.assertFalse(out.exists())
            self.assertEqual(list(root.glob('.eic-stage-*')), [])

    def test_plot_and_eic_cleanup_is_scoped_to_known_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plots = root / 'plots'
            (plots / 'tic').mkdir(parents=True)
            (plots / 'tic' / 'old.png').write_text('old', encoding='utf-8')
            (plots / 'custom').mkdir()
            (plots / 'custom' / 'user.txt').write_text('keep', encoding='utf-8')
            stage = begin_output(plots)
            (stage / 'tic' / 'new.png').parent.mkdir(parents=True)
            (stage / 'tic' / 'new.png').write_text('new', encoding='utf-8')
            commit_output(stage, plots, group='plot')
            self.assertFalse((plots / 'tic' / 'old.png').exists())
            self.assertEqual((plots / 'custom' / 'user.txt').read_text(encoding='utf-8'), 'keep')

            eic_dir = root / 'eic'
            (eic_dir / 'plots').mkdir(parents=True)
            (eic_dir / 'old.csv').write_text('old', encoding='utf-8')
            (eic_dir / 'plots' / 'old.png').write_text('old', encoding='utf-8')
            (eic_dir / 'plots' / 'user.txt').write_text('keep', encoding='utf-8')
            (eic_dir / 'user.csv').write_text('keep', encoding='utf-8')
            (eic_dir / 'peaks.csv').write_text('run\nold\n', encoding='utf-8')
            (eic_dir / 'empty-user-dir').mkdir()
            stage = begin_output(eic_dir)
            (stage / 'new.csv').write_text('new', encoding='utf-8')
            commit_output(stage, eic_dir, group='eic')
            self.assertFalse((eic_dir / 'old.csv').exists())
            self.assertFalse((eic_dir / 'plots' / 'old.png').exists())
            self.assertEqual((eic_dir / 'plots' / 'user.txt').read_text(encoding='utf-8'), 'keep')
            self.assertEqual((eic_dir / 'user.csv').read_text(encoding='utf-8'), 'keep')
            self.assertTrue((eic_dir / 'empty-user-dir').is_dir())


if __name__ == '__main__':
    unittest.main()
