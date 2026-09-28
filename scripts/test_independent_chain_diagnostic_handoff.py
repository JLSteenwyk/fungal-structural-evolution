#!/usr/bin/env python3
"""Synthetic receipt fixtures for four-chain diagnostic provenance checks."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from ancestral_chain_attempt import sha, write_json
from advance_independent_chain_diagnostics import quartet_manifest


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.jobs = []
        self.mapping_hash = 'synthetic-mapping-hash'
        for i in range(4):
            chain = dict(chain_id='synthetic-chain-%d' % i, seed=i+1)
            config = dict(command=['synthetic-only'], seed=i+1, model_input_identity='synthetic-input')
            self.jobs.append(dict(chain=chain, config=config))
            folder = self.root / chain['chain_id']
            attempt = folder / 'attempt-0001'
            attempt.mkdir(parents=True)
            write_json(folder / 'configuration.json', config)
            write_json(attempt / 'command.json', config['command'])
            headers = ['iter', 'scale', 'scale1', 'scale*|T|', 'scale1*|T|', '|T|',
                       'prior', 'likelihood', 'posterior', 'ASRV.Gamma:alpha', 'RS07:rate', 'RS07:meanLength']
            headers += ['F:pi[%s]' % aa for aa in 'ACDEFGHIKLMNPQRSTVWY']
            log = attempt / 'C1.log'
            log.write_text('\t'.join(headers) + '\n')
            receipt = attempt / 'receipt.json'
            digest = hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
            write_json(receipt, dict(status='exited_zero_pending_scientific_validation', exit_code=0,
                configuration_sha256=digest, artifacts={p.name: sha(p) for p in attempt.iterdir()}))
            write_json(folder / 'attempt-0001-sample-audit.json', dict(
                status='all_saved_alignments_and_candidate_nodes_checked', iterations=1000,
                mapping_sha256=self.mapping_hash, attempt_receipt=str(receipt),
                attempt_receipt_sha256=sha(receipt), scalar_log=str(log), scalar_log_sha256=sha(log)))

    def tearDown(self):
        self.temporary.cleanup()

    def test_ready_quartet(self):
        manifest, evidence = quartet_manifest(self.jobs, self.root, 1000, self.mapping_hash)
        self.assertEqual(len(manifest['chains']), 4)
        self.assertEqual(len(manifest['variables']), 26)
        self.assertEqual(len(evidence), 8)
        self.assertEqual(manifest['expected_iterations'], list(range(1001)))

    def test_missing_member_waits(self):
        (self.root / 'synthetic-chain-3/attempt-0001-sample-audit.json').unlink()
        self.assertIsNone(quartet_manifest(self.jobs, self.root, 1000, self.mapping_hash))

    def test_artifact_tampering_rejected(self):
        (self.root / 'synthetic-chain-0/attempt-0001/C1.log').write_text('changed')
        with self.assertRaises(AssertionError):
            quartet_manifest(self.jobs, self.root, 1000, self.mapping_hash)

    def test_input_or_seed_changes_rejected(self):
        self.jobs[0]['config']['seed'] = 100
        with self.assertRaises(AssertionError):
            quartet_manifest(self.jobs, self.root, 1000, self.mapping_hash)


if __name__ == '__main__':
    unittest.main(verbosity=2)
