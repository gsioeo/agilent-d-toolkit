"""Configuration loading, validation and the concentration series."""
import json
import copy
import tempfile
import unittest
from pathlib import Path

import context

from mrm_quant import config

HEADER = ','.join(config.BATCH_COLUMNS)
ANALYTE = {
    'version': 'v1',
    'analytes': [{
        'analyte_id': 'A', 'expected_rt_min': 14.43, 'rt_tolerance_min': 0.08,
        'quantifier': {'precursor_mz': 204., 'product_mz': 93., 'collision_energy_ev': 30.,
                       'polarity': 0, 'time_segment_id': 1, 'scan_method_id': 1},
        'calibration': {'weighting': '1/x', 'intercept': 'free',
                        'concentration_unit': 'mg/mL',
                        'series': {'top_concentration': 1., 'dilution_step': 2.,
                                   'levels': 10}}}]}


def write(directory, name, text):
    path = Path(directory) / name
    path.write_text(text, encoding='utf-8')
    return path


def batch_text(rows):
    return HEADER + '\n' + '\n'.join(rows) + '\n'


def row(run_id='CAL_1', role='calibration', level='S1', concentration='',
        unit='', dilution='1', include='true', reason='', analyte='A',
        dataset=None, batch='B1'):
    return ','.join([run_id, dataset or ('%s.d' % run_id), batch, analyte, role, level,
                     concentration, unit, '', dilution, '', '', include, reason])


class SeriesTests(unittest.TestCase):
    def test_two_fold_series_from_one_milligram_per_millilitre(self):
        series = config.concentration_series(top_concentration=1., dilution_step=2.,
                                             levels=10)
        self.assertEqual(series['S1'], 1.)
        self.assertEqual(series['S10'], 1. / 512)
        self.assertEqual(len(series), 10)

    def test_series_parameters_are_validated(self):
        cases = [dict(top_concentration=0., dilution_step=2., levels=10),
                 dict(top_concentration=1., dilution_step=1., levels=10),
                 dict(top_concentration=1., dilution_step=2., levels=0)]
        for arguments in cases:
            with self.subTest(**arguments), self.assertRaisesRegex(ValueError, 'invalid_series'):
                config.concentration_series(**arguments)

    def test_level_identifier_format_is_configurable(self):
        series = config.concentration_series(top_concentration=100., dilution_step=10.,
                                            levels=3, level_id_format='L%02d')
        self.assertEqual(sorted(series), ['L01', 'L02', 'L03'])
        self.assertAlmostEqual(series['L03'], 1.)


class BatchTests(unittest.TestCase):
    def load(self, rows, **options):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(tmp, 'batch.csv', batch_text(rows))
            options.setdefault('require_datasets', False)
            return config.load_batch(path, **options)

    def test_unknown_column_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(tmp, 'batch.csv', 'run_id,dataset_path,batch_id,analyte_id,role,extra\n')
            with self.assertRaisesRegex(ValueError, 'invalid_batch_columns'):
                config.load_batch(path, require_datasets=False)

    def test_missing_required_column_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(tmp, 'batch.csv', 'run_id,dataset_path,batch_id,analyte_id\n')
            with self.assertRaisesRegex(ValueError, 'invalid_batch_columns'):
                config.load_batch(path, require_datasets=False)

    def test_duplicate_run_and_analyte_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'duplicate_run_analyte'):
            self.load([row(), row()])

    def test_role_must_be_explicit(self):
        with self.assertRaisesRegex(ValueError, 'invalid_role'):
            self.load([row(role='')])
        with self.assertRaisesRegex(ValueError, 'invalid_role'):
            self.load([row(role='standard')])

    def test_exclusion_needs_a_reason(self):
        with self.assertRaisesRegex(ValueError, 'exclusion_reason_required'):
            self.load([row(include='false')])
        batch = self.load([row(include='false', reason='bad injection')])
        self.assertFalse(batch['rows'][0]['include'])

    def test_dilution_factor_defaults_to_one_and_must_be_positive(self):
        batch = self.load([row(dilution='')])
        self.assertEqual(batch['rows'][0]['dilution_factor'], 1.)
        with self.assertRaisesRegex(ValueError, 'invalid_dilution'):
            self.load([row(dilution='0')])

    def test_same_run_with_two_dataset_paths_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'inconsistent_dataset_path'):
            self.load([row(analyte='A', dataset='one.d'), row(analyte='B', dataset='two.d')])

    def test_missing_dataset_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(tmp, 'batch.csv', batch_text([row()]))
            with self.assertRaisesRegex(ValueError, 'dataset_not_found'):
                config.load_batch(path)


class AnalyteAndCrossCheckTests(unittest.TestCase):
    def configuration(self, rows, analytes=None):
        with tempfile.TemporaryDirectory() as tmp:
            batch_path = write(tmp, 'batch.csv', batch_text(rows))
            analytes_path = write(tmp, 'analytes.json',
                                  json.dumps(analytes or ANALYTE, ensure_ascii=False))
            return config.load_configuration(batch_path, analytes_path,
                                             require_datasets=False)

    def test_series_fills_missing_concentrations(self):
        rows = [row('CAL_1', level='S1'), row('CAL_2', level='S2'), row('CAL_3', level='S3')]
        loaded = self.configuration(rows)
        values = [item['concentration'] for item in loaded['batch']['rows']]
        self.assertEqual(values, [1., 0.5, 0.25])
        self.assertEqual({item['concentration_source'] for item in loaded['batch']['rows']},
                         {'series'})

    def test_internal_standard_must_reference_external_analyte(self):
        analytes = copy.deepcopy(ANALYTE)
        analytes['analytes'][0]['response'] = {
            'mode': 'internal', 'internal_standard_id': 'IS'}
        is_entry = copy.deepcopy(analytes['analytes'][0])
        is_entry['analyte_id'] = 'IS'
        is_entry['response'] = {'mode': 'external', 'internal_standard_id': None}
        analytes['analytes'].append(is_entry)
        loaded = self.configuration([row('CAL_1')], analytes=analytes)
        self.assertEqual(loaded['analytes']['analytes']['A']['response']['mode'], 'internal')
        bad = copy.deepcopy(analytes)
        bad['analytes'][0]['response']['internal_standard_id'] = 'MISSING'
        with self.assertRaisesRegex(ValueError, 'unknown_internal_standard'):
            self.configuration([row('CAL_1')], analytes=bad)
        self.assertEqual({item['concentration_unit'] for item in loaded['batch']['rows']},
                         {'mg/mL'})

    def test_written_concentration_must_agree_with_the_series(self):
        rows = [row('CAL_1', level='S1', concentration='1', unit='mg/mL'),
                row('CAL_2', level='S2', concentration='0.4', unit='mg/mL'),
                row('CAL_3', level='S3', concentration='0.25', unit='mg/mL')]
        with self.assertRaisesRegex(ValueError, 'concentration_conflict'):
            self.configuration(rows)

    def test_matching_written_concentration_is_kept_with_its_source(self):
        rows = [row('CAL_1', level='S1', concentration='1', unit='mg/mL'),
                row('CAL_2', level='S2', concentration='0.5', unit='mg/mL'),
                row('CAL_3', level='S3', concentration='0.25', unit='mg/mL')]
        loaded = self.configuration(rows)
        self.assertEqual({item['concentration_source'] for item in loaded['batch']['rows']},
                         {'batch_table'})

    def test_unit_must_match_the_analyte(self):
        rows = [row('CAL_1', level='S1', unit='ng/mL')]
        with self.assertRaisesRegex(ValueError, 'unit_mismatch'):
            self.configuration(rows)

    def test_unknown_level_and_missing_concentration_are_refused(self):
        with self.assertRaisesRegex(ValueError, 'unknown_level'):
            self.configuration([row('CAL_1', level='S99')])
        analytes = json.loads(json.dumps(ANALYTE))
        analytes['analytes'][0]['calibration']['series'] = None
        with self.assertRaisesRegex(ValueError, 'missing_concentration'):
            self.configuration([row('CAL_1', level='S1')], analytes)

    def test_unknown_analyte_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'unknown_analyte'):
            self.configuration([row('CAL_1', analyte='Z')])

    def test_analyte_entry_is_validated(self):
        analytes = json.loads(json.dumps(ANALYTE))
        del analytes['analytes'][0]['expected_rt_min']
        with self.assertRaisesRegex(ValueError, 'invalid_analytes'):
            self.configuration([row()], analytes)
        analytes = json.loads(json.dumps(ANALYTE))
        analytes['analytes'][0]['integration'] = {'baseline': 'none', 'typo': 1}
        with self.assertRaisesRegex(ValueError, 'invalid_analytes'):
            self.configuration([row()], analytes)
        analytes = json.loads(json.dumps(ANALYTE))
        del analytes['analytes'][0]['quantifier']['product_mz']
        with self.assertRaisesRegex(ValueError, 'invalid_transition'):
            self.configuration([row()], analytes)

    def test_configuration_hash_covers_both_files(self):
        rows = [row('CAL_1', level='S1'), row('CAL_2', level='S2'), row('CAL_3', level='S3')]
        first = self.configuration(rows)
        second = self.configuration(rows[:2] + [row('CAL_3', level='S4')])
        self.assertNotEqual(first['config_hash'], second['config_hash'])
        self.assertEqual(len(first['config_hash']), 64)

    def test_shipped_templates_load(self):
        templates = context.REPO / 'mrm_quant' / 'templates'
        analytes = config.load_analytes(templates / 'analytes.example.json')
        self.assertIn('<analyte-id>', analytes['analytes'])
        batch = config.load_batch(templates / 'batch.example.csv', require_datasets=False)
        self.assertEqual({item['role'] for item in batch['rows']},
                         {'calibration', 'blank', 'qc', 'unknown'})


if __name__ == '__main__':
    unittest.main()
