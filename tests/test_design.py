import unittest
from helixbusters.design import validate_design


class TestExperimentalDesign(unittest.TestCase):
    def rows(self):
        return [{'sample': f'{d}_{g}', 'group': g, 'replicate': f'R{i}', 'donor': d}
                for g in ('ACUTE', 'CHRONIC') for i, d in enumerate(('HD1', 'HD2', 'HD3'), 1)]

    def test_complete_pairing_uses_donor_not_replicate_label(self):
        rows = self.rows()
        rows[-1]['replicate'] = 'different_label'
        report = validate_design(rows, 'paired')
        self.assertEqual(report['declared_donors'], 3)
        self.assertEqual(report['model_formula'], '~ donor + condition')
        self.assertEqual(report['sample_metadata'][-1]['donor'], 'HD3')

    def test_paired_requires_explicit_and_complete_donors(self):
        rows = self.rows()
        rows[0]['donor'] = ''
        with self.assertRaisesRegex(ValueError, 'requires donor'):
            validate_design(rows, 'paired')
        with self.assertRaisesRegex(ValueError, 'Incomplete pairing'):
            validate_design(self.rows()[:-1], 'paired')

    def test_unpaired_rejects_reused_donors_and_accepts_distinct_donors(self):
        with self.assertRaisesRegex(ValueError, 'reuses donors'):
            validate_design(self.rows(), 'unpaired')
        rows = self.rows()
        for row in rows:
            row['donor'] = row['sample']
        self.assertEqual(validate_design(rows, 'unpaired')['model_formula'], '~ condition')

    def test_legacy_replicate_labels_do_not_infer_pairing(self):
        rows = self.rows()
        for row in rows:
            del row['donor']
        self.assertEqual(validate_design(rows)['design'], 'unspecified')
        self.assertIsNone(validate_design(rows)['model_formula'])
        self.assertEqual(validate_design(rows, 'unpaired')['declared_donors'], 0)

    def test_duplicate_donor_in_condition_is_rejected(self):
        rows = self.rows()
        rows[1]['donor'] = 'HD1'
        with self.assertRaisesRegex(ValueError, 'Duplicate donor'):
            validate_design(rows, 'paired')
