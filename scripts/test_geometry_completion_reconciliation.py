#!/usr/bin/env python3
"""Check full-set reconciliation and legitimate longer rotation degeneracy."""
import unittest
from record_primary_diagnostic_geometry_completion_v2 import reconcile, launch_plan


def row(identity, length, status='degenerate_at_numeric_tolerance', rmsd='within_printed_rounding'):
    return dict(pair_key=identity, mask='full', order=0, aligned_length=length,
        geometry_status=status, rmsd_status=rmsd)


class ReconciliationTests(unittest.TestCase):
    def test_legacy_launch_binding(self):
        launch = dict(cmdline=['python', 'script.py', '--plan', 'old.json'])
        self.assertEqual(launch_plan(launch), 'old.json')
        with self.assertRaises(ValueError):launch_plan(dict(launch, plan='wrong.json'))

    def setUp(self):
        self.short = row('short', 2, rmsd='outside_printed_rounding')
        self.long = row('collinear', 20)
        self.unique = row('unique', 30, 'unique_at_numeric_tolerance')
        self.all = [self.short, self.long, self.unique]

    def run_rows(self, geometry=None, numeric=None, short=None, audit=None):
        return reconcile(self.all if geometry is None else geometry,
            self.all if numeric is None else numeric,
            [self.short] if short is None else short,
            [self.short, self.long] if audit is None else audit)

    def test_long_degeneracy_retained(self):
        result = self.run_rows()
        self.assertEqual(result['degenerate_long_alignments'], 1)
        self.assertEqual(result['longer_degenerate_mappings'][0]['pair_key'], 'collinear')
        self.assertEqual(result['numerically_unique_rotations'], 1)
        self.assertEqual(result['rmsd_classification_counts']['outside_printed_rounding'], 1)

    def test_missing_row_rejected(self):
        with self.assertRaises(ValueError):self.run_rows(geometry=self.all[:2])

    def test_duplicate_row_rejected(self):
        with self.assertRaises(ValueError):self.run_rows(geometry=self.all+[self.long])

    def test_changed_rmsd_exclusion_rejected(self):
        changed = dict(self.short, rmsd_status='within_printed_rounding')
        with self.assertRaises(ValueError):self.run_rows(geometry=[changed, self.long, self.unique])

    def test_short_unique_rotation_rejected(self):
        changed = dict(self.short, geometry_status='unique_at_numeric_tolerance')
        with self.assertRaises(ValueError):self.run_rows(geometry=[changed, self.long, self.unique])

    def test_missing_analytic_census_rejected(self):
        with self.assertRaises(ValueError):self.run_rows(short=[])

    def test_missing_long_audit_rejected(self):
        with self.assertRaises(ValueError):self.run_rows(audit=[self.short])

    def test_duplicate_audit_rejected(self):
        with self.assertRaises(ValueError):self.run_rows(audit=[self.short, self.long, self.long])


if __name__ == '__main__':unittest.main()
