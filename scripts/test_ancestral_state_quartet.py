#!/usr/bin/env python3
import copy
import json
from pathlib import Path
import unittest
import numpy as np
from ancestral_chain_attempt import sha,write_json
import test_independent_chain_diagnostic_handoff as handoff_tests
from read_ancestral_state_quartet import read_quartet,validate_arrays,ALPHABET


class StateQuartetTests(handoff_tests.HandoffTests):
    def setUp(self):
        super().setUp();self.state_root=self.root/'state-extraction'
        alignment=self.root/'observed.fasta';alignment.write_text('>a\nAA\n>b\nAA\n')
        self.coordinates=dict(alphabet=ALPHABET,nodes=['n0','n1','n2','n3'],tips=[dict(tip='a',length=2),dict(tip='b',length=2)],input_alignment=str(alignment),input_alignment_sha256=sha(alignment),ordering='synthetic fixture coordinates')
        self.arrays=[]
        for i,job in enumerate(self.jobs):
            chain=job['chain'];chain.update(alignment=str(alignment),alignment_sha256=sha(alignment))
            ap=self.root/chain['chain_id']/'attempt-0001-sample-audit.json';audit=json.loads(ap.read_text())
            audit['candidate_samples']=[dict(iteration=t,source_node=n) for t in range(0,1001,10) for n in self.coordinates['nodes']];write_json(ap,audit)
            folder=self.state_root/chain['chain_id'];attempt=folder/'processing/attempt-0001';out=attempt/'states';target=out/chain['chain_id'];target.mkdir(parents=True)
            values=np.random.default_rng(i).integers(0,3,size=(101,4,4),dtype=np.uint8);iterations=np.arange(0,1001,10)
            saved=dict(states=values,iterations=iterations,unanchored_residue_counts=np.zeros((101,4),dtype=np.int32))
            for cutoff in [250,500]:saved[f'counts_after_{cutoff}']=np.stack([np.count_nonzero(values[iterations>cutoff]==s,axis=0) for s in range(len(ALPHABET))],axis=-1).astype(np.uint16)
            self.arrays.append(saved);np.savez_compressed(target/'states.npz',**saved);write_json(target/'coordinates.json',self.coordinates)
            summary=dict(chain=chain['chain_id'],samples=101,state_observations=values.size,source_audit=str(ap),source_audit_sha256=sha(ap),artifacts={str(p):sha(p) for p in target.iterdir()})
            rp=out/'receipt.json';write_json(rp,dict(status='complete_independently_checked_anchored_state_traces',chains=1,pins={},summaries=[summary]))
            ar=attempt/'receipt.json';write_json(ar,dict(exit_code=0,artifacts={str(p.relative_to(attempt)):sha(p) for p in attempt.rglob('*') if p.is_file()}))
            write_json(folder/'disposition.json',dict(status='state_trace_complete_not_posterior_qualification',receipt=str(rp),receipt_sha256=sha(rp),attempt_receipt=str(ar),attempt_receipt_sha256=sha(ar)))

    def read(self):return read_quartet(self.jobs,self.root,self.state_root,1000,self.mapping_hash)

    def test_complete_arrays_bind_to_quartet(self):
        result=self.read();self.assertEqual(result['values'].shape,(4,101,4,4));self.assertEqual(len(result['chains']),4)
        for i,saved in enumerate(self.arrays):np.testing.assert_array_equal(result['values'][i],saved['states'])

    def test_missing_extraction_waits(self):
        (self.state_root/'synthetic-chain-3/disposition.json').unlink();self.assertIsNone(self.read())

    def test_state_file_tampering_rejected(self):
        p=next(self.state_root.rglob('states.npz'));p.write_bytes(b'altered')
        with self.assertRaises(AssertionError):self.read()

    def test_counts_schedule_and_nodes_rejected(self):
        saved=copy.deepcopy(self.arrays[0]);saved['counts_after_250'][0,0,0]+=1
        with self.assertRaises(AssertionError):validate_arrays(saved,self.coordinates,list(range(0,1001,10)))
        saved=copy.deepcopy(self.arrays[0]);saved['iterations'][1]=11
        with self.assertRaises(AssertionError):validate_arrays(saved,self.coordinates,list(range(0,1001,10)))
        coords=copy.deepcopy(self.coordinates);coords['nodes'][1]=coords['nodes'][0]
        with self.assertRaises(AssertionError):validate_arrays(self.arrays[0],coords,list(range(0,1001,10)))


if __name__=='__main__':unittest.main()
