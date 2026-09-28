#!/usr/bin/env python3
"""Known mixing failures and passing synthetic reference distributions."""
import unittest
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from ancestral_chain_diagnostics import diagnose


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(20260927)

    def test_independent_normal(self):
        result = diagnose(self.rng.normal(size=(4, 4000)))
        self.assertEqual(result['status'], 'passes_scalar_screen_only')

    def test_shifted_chain(self):
        draws = self.rng.normal(size=(4, 4000))
        draws[0] += 3
        self.assertEqual(diagnose(draws)['status'], 'scalar_mixing_requires_review')

    def test_scale_mismatch(self):
        draws = self.rng.normal(size=(4, 4000))
        draws[0] *= 10
        self.assertGreater(diagnose(draws)['rhat'], 1.01)

    def test_slow_autocorrelation(self):
        draws = self.rng.normal(size=(4, 4000))
        for i in range(1, draws.shape[1]):
            draws[:, i] += .999 * draws[:, i - 1]
        result = diagnose(draws)
        self.assertEqual(result['status'], 'scalar_mixing_requires_review')
        self.assertLess(result['bulk_ess'], 400)

    def test_degenerate_inputs(self):
        self.assertEqual(diagnose(np.ones((4, 100)))['status'], 'constant_chain_requires_review')
        self.assertEqual(diagnose(np.ones((4, 3)))['status'], 'insufficient_retained_draws')
        self.assertEqual(diagnose(np.full((4, 100), np.nan))['status'], 'nonfinite_draws')
        with self.assertRaises(ValueError):
            diagnose(np.ones((3, 100)))

    def test_log_manifest_integration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chains = []
            for chain in range(4):
                path = root / ('chain%d.tsv' % chain)
                with path.open('w') as handle:
                    writer = csv.writer(handle, delimiter='\t')
                    writer.writerow(['iter', 'normal', 'shifted'])
                    for i, value in enumerate(self.rng.normal(size=2001)):
                        writer.writerow([i, value, value + (3 if chain == 0 else 0)])
                chains.append(dict(chain_id=str(chain), seed=chain + 1,
                    model_input_identity='synthetic-test-only', log=str(path),
                    log_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            manifest = dict(chains=chains, expected_iterations=list(range(2001)),
                            discard_through_iteration=500, variables=['normal', 'shifted'])
            path = root / 'manifest.json'
            path.write_text(json.dumps(manifest))
            command = [sys.executable, str(Path(__file__).with_name('ancestral_chain_diagnostics.py')),
                       '--manifest', str(path), '--output', str(root / 'result')]
            subprocess.run(command, check=True, capture_output=True)
            result = json.loads((root / 'result/diagnostics.json').read_text())
            self.assertEqual(result['variables']['normal']['draws_per_chain'], 1500)
            self.assertEqual(result['variables']['normal']['status'], 'passes_scalar_screen_only')
            self.assertEqual(result['variables']['shifted']['status'], 'scalar_mixing_requires_review')
            manifest['chains'][1]['seed'] = manifest['chains'][0]['seed']
            path.write_text(json.dumps(manifest))
            command[-1] = str(root / 'invalid-seeds')
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertFalse((root / 'invalid-seeds').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
