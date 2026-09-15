"""The three commands end to end, including output-directory protection."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

import test_reader_adapter as adapter

from mrm_quant import __main__ as cli


def read_csv(path):
    with open(path, encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


class OutputProtectionTests(unittest.TestCase):
    def test_existing_directory_is_refused_and_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset = adapter.synthetic_dataset(tmp)
            existing = Path(tmp) / 'results'
            existing.mkdir()
            code = cli.main(['inspect', '--source', str(dataset), '--out', str(existing)])
            self.assertEqual(code, 1)
            self.assertTrue((existing / 'failure.json').is_file())

    def test_output_inside_a_dataset_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset = adapter.synthetic_dataset(tmp)
            code = cli.main(['inspect', '--source', str(dataset),
                             '--out', str(dataset / 'results')])
            self.assertEqual(code, 1)
            self.assertFalse((dataset / 'results').exists())


class SyntheticInspectTests(unittest.TestCase):
    def test_inspect_writes_tables_for_a_dataset_without_a_method(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset = adapter.synthetic_dataset(tmp)
            out = Path(tmp) / 'inspection'
            self.assertEqual(cli.main(['inspect', '--source', str(dataset),
                                       '--out', str(out)]), 0)
            runs = read_csv(out / 'runs.csv')
            self.assertEqual(len(runs), 1)
            self.assertEqual(runs[0]['intensity_unit'], 'stored_intensity')
            self.assertEqual(read_csv(out / 'channels.csv'), [])
            provenance = json.loads((out / 'provenance.json').read_text(encoding='utf-8'))
            self.assertEqual(provenance['command'], 'inspect')
            self.assertEqual(provenance['area_unit'], 'stored_intensity*min')


if __name__ == '__main__':
    unittest.main()
