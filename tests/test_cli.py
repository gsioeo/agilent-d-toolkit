"""The three commands end to end, including output-directory protection."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

import context
import september
import test_reader_adapter as adapter

from mrm_quant import __main__ as cli
from mrm_quant import config

CALIBRATION_RUNS = tuple('STD_S%d' % level for level in range(1, 11))
SAMPLE_RUNS = ('S21',)


def read_csv(path):
    with open(path, encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def batch_text(data_root):
    lines = [','.join(config.BATCH_COLUMNS)]
    for index, name in enumerate(CALIBRATION_RUNS, start=1):
        lines.append(','.join([name, 'TSJ-0907/LQ/%s.d' % name, september.BATCH_ID,
                               september.ANALYTE_ID, 'calibration', 'S%d' % index,
                               '', '', 'vial concentration', '1', '', '', 'true', '']))
    for name in SAMPLE_RUNS:
        lines.append(','.join([name, 'TSJ-0907/LQ/%s.d' % name, september.BATCH_ID,
                               september.ANALYTE_ID, 'unknown', '', '', '',
                               '', '1', '', '', 'true', '']))
    return '\n'.join(lines) + '\n'


def analytes_document():
    def transition(product, label):
        return {'precursor_mz': 204., 'product_mz': product, 'collision_energy_ev': 30.,
                'polarity': 0, 'time_segment_id': 1, 'scan_method_id': 1,
                'channel_label': label}
    return {'version': 'v1', 'analytes': [{
        'analyte_id': september.ANALYTE_ID,
        'identity_confirmed': False,
        'expected_rt_min': september.EXPECTED_RT_MIN,
        'rt_tolerance_min': september.RT_TOLERANCE_MIN,
        'search_margin_min': september.SEARCH_MARGIN_MIN,
        'quantifier': transition(september.QUANTIFIER_PRODUCT_MZ, '204to93'),
        'qualifiers': [transition(81., '204to81')],
        'integration': {'baseline': september.BASELINE,
                        'min_relative_height': september.FIND_OPTIONS['min_relative_height'],
                        'min_separation_min': september.FIND_OPTIONS['min_separation_min']},
        'calibration': {'weighting': september.WEIGHTING, 'intercept': september.INTERCEPT,
                        'concentration_unit': september.CONCENTRATION_UNIT,
                        'series': {'top_concentration': september.TOP_CONCENTRATION,
                                   'dilution_step': september.DILUTION_STEP,
                                   'levels': september.LEVELS}},
        'response': {'mode': 'external'},
        'qc': {'blank_run_id': None}}]}


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


class SeptemberEndToEndTests(unittest.TestCase):
    def setUp(self):
        if not september.dataset('STD_S1').is_dir():
            self.skipTest('Local raw dataset unavailable')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.data_root = str(context.DATA_ROOT)
        self.batch = self.root / 'batch.csv'
        self.batch.write_text(batch_text(self.data_root), encoding='utf-8')
        self.analytes = self.root / 'analytes.json'
        self.analytes.write_text(json.dumps(analytes_document(), ensure_ascii=False),
                                 encoding='utf-8')

    def run_command(self, name, out, *extra):
        return cli.main([name, '--batch', str(self.batch), '--analytes', str(self.analytes),
                         '--data-root', self.data_root, '--out', str(out), *extra])

    def test_inspect_ranks_the_strongest_channel_first(self):
        out = self.root / 'inspection'
        self.assertEqual(cli.main(['inspect', '--source',
                                   str(september.dataset('STD_S1')), '--out', str(out),
                                   '--rt-range', '14.2 14.7']), 0)
        rows = read_csv(out / 'channel_response.csv')
        first = [row for row in rows if row['rank_in_run'] == '1']
        self.assertEqual(len(first), 1)
        self.assertEqual(float(first[0]['product_mz']), september.QUANTIFIER_PRODUCT_MZ)
        provenance = json.loads((out / 'provenance.json').read_text(encoding='utf-8'))
        self.assertIn('quantify never chooses a channel by response',
                      provenance['note'])

    def test_extract_writes_one_trace_per_configured_channel(self):
        out = self.root / 'traces'
        self.assertEqual(self.run_command('extract', out, '--no-checksums'), 0)
        traces = read_csv(out / 'traces.csv')
        self.assertEqual(len(traces), 2 * (len(CALIBRATION_RUNS) + len(SAMPLE_RUNS)))
        self.assertEqual({row['points'] for row in traces}, {'2454'})
        first = read_csv(out / 'traces' / 'STD_S1' / 'Ses1_204to93.csv')
        self.assertEqual(len(first), 2454)
        self.assertEqual(float(first[0]['intensity']), 1121.11181640625)
        self.assertEqual(first[0]['point_status'], 'ok')

    def test_quantify_reports_concentration_equation_and_r2(self):
        out = self.root / 'quant'
        self.assertEqual(self.run_command('quantify', out, '--no-checksums', '--no-plots'), 0)
        results = {row['run_id']: row for row in read_csv(out / 'results.csv')}
        self.assertEqual(results['S21']['status'], 'ok')
        self.assertGreater(float(results['S21']['reported_concentration']), 0.)
        self.assertEqual(results['S21']['concentration_unit'], 'mg/mL')
        self.assertEqual(results['S21']['validated'], 'False')
        self.assertTrue(results['S21']['calibration_equation'].startswith('A = '))
        self.assertGreater(float(results['S21']['r2']), 0.9)

        models = json.loads((out / 'calibration_models.json').read_text(encoding='utf-8'))
        model = models['models'][0]['model']
        self.assertEqual(model['n_levels'], len(CALIBRATION_RUNS))
        self.assertEqual(model['weighting'], september.WEIGHTING)
        self.assertEqual(model['concentration_unit'], 'mg/mL')
        self.assertIn('Descriptive only', model['r2_definition'])

        points = read_csv(out / 'calibration_points.csv')
        self.assertEqual([point['level_id'] for point in points],
                         ['S%d' % index for index in range(1, len(CALIBRATION_RUNS) + 1)])
        self.assertEqual({point['concentration_source'] for point in points}, {'series'})

        qc_rows = {row['run_id']: row for row in read_csv(out / 'qc.csv')}
        self.assertEqual(qc_rows['S21']['blank_status'], 'not_evaluated')
        self.assertEqual(qc_rows['S21']['independent_qc_status'], 'not_evaluated')
        self.assertIn('204to81=not_evaluated', qc_rows['S21']['qualifier_summary'])

        provenance = json.loads((out / 'provenance.json').read_text(encoding='utf-8'))
        self.assertEqual(provenance['configuration']['batch_ids'], [september.BATCH_ID])
        self.assertEqual(len(provenance['configuration']['config_hash']), 64)
        self.assertEqual(provenance['analyte_settings']['Ses1']['expected_rt_min'],
                         september.EXPECTED_RT_MIN)
        self.assertFalse(provenance['analyte_settings']['Ses1']['identity_confirmed'])

    def test_quantify_refuses_a_batch_with_mixed_methods(self):
        mixed = self.batch.read_text(encoding='utf-8').replace(
            'TSJ-0907/LQ/S21.d', '20260819-tsj/TSJ-0819-1.d')
        self.batch.write_text(mixed, encoding='utf-8')
        if not context.dataset_path('20260819-tsj/TSJ-0819-1.d').is_dir():
            self.skipTest('August dataset unavailable')
        out = self.root / 'mixed'
        self.assertEqual(self.run_command('quantify', out, '--no-checksums', '--no-plots'), 1)
        failure = json.loads((out / 'failure.json').read_text(encoding='utf-8'))
        self.assertIn('method_mismatch', failure['error'])


if __name__ == '__main__':
    unittest.main()
