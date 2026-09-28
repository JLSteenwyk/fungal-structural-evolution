#!/usr/bin/env python3
import gzip
import json
from pathlib import Path
import subprocess
import unittest
import numpy as np
import test_ancestral_state_quartet as fixtures
from ancestral_chain_attempt import sha,write_json


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.StateQuartetTests();self.fixture.setUp()
        self.root=self.fixture.root;self.quartet=self.fixture.read();q=self.quartet
        self.arrays=self.root/'quartet.npz'
        np.savez_compressed(self.arrays,values=q['values'],iterations=q['iterations'],unanchored_residue_counts=np.stack([a['unanchored_residue_counts'] for a in self.fixture.arrays]))
        self.manifest=self.root/'manifest.json'
        self.payload=dict(status='provenance_checked_state_quartet',chains=q['chains'],coordinates=q['coordinates'],expected_iterations=q['iterations'].tolist(),arrays=str(self.arrays),arrays_sha256=sha(self.arrays),evidence=q['evidence'])
        write_json(self.manifest,self.payload)
        self.command=['SOFTWARE/ancestral-diagnostics-20260927/bin/python','scripts/report_ancestral_state_quartet.py','--manifest',str(self.manifest),'--output',str(self.root/'reports')]

    def tearDown(self):self.fixture.tearDown()

    def test_full_reader_to_report(self):
        completed=subprocess.run(self.command,capture_output=True,text=True)
        self.assertEqual(completed.returncode,0,completed.stderr)
        root=self.root/'reports';receipt=json.loads((root/'receipt.json').read_text())
        self.assertEqual(set(receipt['outputs']),{'250','500'})
        for cutoff,draws in [(250,75),(500,50)]:
            folder=root/f'discard-{cutoff}';summary=json.loads((folder/'summary.json').read_text())
            self.assertEqual(summary['coordinates'],16);self.assertEqual(summary['retained_samples_per_chain'],draws)
            self.assertEqual(sum(summary['coordinate_status_counts'].values()),16)
            with np.load(folder/'patterns.npz',allow_pickle=False) as saved:
                patterns=saved['patterns'];mapping=saved['coordinate_pattern_ids']
            with gzip.open(folder/'diagnostics.jsonl.gz','rt') as handle:reports=[json.loads(line) for line in handle]
            self.assertEqual(sum(r['coordinate_multiplicity'] for r in reports),16)
            for n in range(4):
                for anchor in range(4):
                    pattern=patterns[mapping[n,anchor]]
                    np.testing.assert_array_equal(pattern,self.quartet['values'][:,self.quartet['iterations']>cutoff,n,anchor])
            self.assertTrue(all(r['diagnostic']['indicators']['-']['screen']['status']=='state_not_observed_probability_unresolved' for r in reports))
        for name,h in receipt['artifacts'].items():self.assertEqual(sha(root/name),h)

    def test_duplicate_seed_rejected_before_output(self):
        self.payload['chains'][1]['seed']=self.payload['chains'][0]['seed'];write_json(self.manifest,self.payload)
        self.assertNotEqual(subprocess.run(self.command,capture_output=True).returncode,0)
        self.assertFalse((self.root/'reports').exists())

    def test_source_array_tampering_rejected_before_output(self):
        self.arrays.write_bytes(b'altered')
        self.assertNotEqual(subprocess.run(self.command,capture_output=True).returncode,0)
        self.assertFalse((self.root/'reports').exists())


if __name__=='__main__':unittest.main()
