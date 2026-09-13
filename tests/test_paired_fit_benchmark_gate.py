"""Ensure failed/changed prerequisites cannot authorize downstream execution."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('gate', Path(__file__).resolve().parents[1] / 'scripts/advance_paired_fit_benchmark.py')
gate = importlib.util.module_from_spec(spec); spec.loader.exec_module(gate)


class Gates(unittest.TestCase):
    def test_changed_pin_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); p = root / 'input'; p.write_text('original')
            c = {'pinned_files': {'input': gate.sha(p)}}
            p.write_text('changed')
            with self.assertRaises(ValueError): gate.verify_pins(c, root)

    def test_missing_or_duplicate_completion_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); fit = root / 'fit'; fit.mkdir(); (fit / 'config.json').write_text('{}')
            c = {'fits': 'fit', 'expected_markers': 2}
            with self.assertRaises(FileNotFoundError): gate.completed_fit(c, root)
            (fit / 'receipt.json').write_text(json.dumps({'status': 'complete_matched_topology_point_estimates', 'config_sha256': gate.sha(fit / 'config.json'), 'results': [{'marker': 'a'}, {'marker': 'a'}]}))
            with self.assertRaises(ValueError): gate.completed_fit(c, root)

    def test_geometry_from_other_input_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for folder in ['audit', 'geometry', 'inputs']:
                (root / folder).mkdir(); (root / folder / 'receipt.json').write_text('{}')
            (root / 'audit/receipt.json').write_text(json.dumps({'status': 'passed_complete_paired_grid_character_and_sampled_geometry_readback', 'source_receipts': {'comparisons': gate.sha(root / 'geometry/receipt.json'), 'inputs': 'wrong'}}))
            with self.assertRaises(ValueError): gate.geometry_gate({'geometry_audit': 'audit', 'geometry': 'geometry', 'inputs': 'inputs'}, root)


if __name__ == '__main__': unittest.main()
