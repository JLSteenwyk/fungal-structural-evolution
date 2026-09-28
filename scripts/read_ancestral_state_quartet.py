#!/usr/bin/env python3
"""Bind state arrays to four checked inference chains before diagnostics."""
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from ancestral_chain_attempt import sha
from advance_independent_chain_diagnostics import quartet_manifest

ALPHABET='ACDEFGHIKLMNPQRSTVWYX-'


def validate_arrays(saved,coordinates,expected):
    values=saved['states'];iterations=saved['iterations'];free=saved['unanchored_residue_counts']
    assert coordinates['alphabet']==ALPHABET
    nodes=coordinates['nodes'];tips=coordinates['tips']
    assert len(nodes)==len(set(nodes))==4 and nodes==sorted(nodes)
    assert tips and len({r['tip'] for r in tips})==len(tips)
    assert [r['tip'] for r in tips]==sorted(r['tip'] for r in tips)
    assert all(isinstance(r['length'],int) and r['length']>0 for r in tips)
    assert np.array_equal(iterations,expected)
    assert values.dtype==np.uint8 and values.shape==(len(expected),4,sum(r['length'] for r in tips))
    assert np.all(values<len(ALPHABET))
    assert free.shape==(len(expected),4) and np.issubdtype(free.dtype,np.integer) and np.all(free>=0)
    for cutoff in [250,500]:
        counts=saved[f'counts_after_{cutoff}'];retained=values[iterations>cutoff]
        assert np.issubdtype(counts.dtype,np.integer) and counts.shape==values.shape[1:]+(len(ALPHABET),)
        assert np.all(counts>=0) and np.all(counts.sum(axis=-1)==len(retained))
        for state in range(len(ALPHABET)):
            assert np.array_equal(counts[...,state],np.count_nonzero(retained==state,axis=0))
    return values.copy()


def read_quartet(jobs,producer_root,state_root,iterations,mapping_hash):
    """Caller must bind jobs/paths and code to verified full producer/extractor plans."""
    assert len(jobs)==len({j['chain']['chain_id'] for j in jobs})==4
    assert len({j['chain']['seed'] for j in jobs})==4
    assert len({j['config']['model_input_identity'] for j in jobs})==1
    expected=list(range(0,iterations+1,10));producer_root=Path(producer_root);state_root=Path(state_root)
    ready=quartet_manifest(jobs,producer_root,iterations,mapping_hash)
    if ready is None:return None
    scalar,evidence=ready
    dispositions=[state_root/j['chain']['chain_id']/'disposition.json' for j in jobs]
    if any(not p.exists() for p in dispositions):return None
    arrays=[];common=None
    for job,dp in zip(jobs,dispositions):
        chain=job['chain'];d=json.loads(dp.read_text())
        if d['status']!='state_trace_complete_not_posterior_qualification':return None
        rp=Path(d['receipt']);assert rp.resolve().is_relative_to(dp.parent.resolve())
        assert sha(rp)==d['receipt_sha256'];report=json.loads(rp.read_text())
        assert report['status']=='complete_independently_checked_anchored_state_traces' and report['chains']==1
        attempt=Path(d['attempt_receipt']);assert sha(attempt)==d['attempt_receipt_sha256']
        ar=json.loads(attempt.read_text());assert ar['exit_code']==0
        for p,h in ar['artifacts'].items():assert sha(attempt.parent/p)==h
        for p,h in report['pins'].items():assert sha(p)==h
        summary=report['summaries'][0];assert summary['chain']==chain['chain_id']
        ap=next((producer_root/chain['chain_id']).glob('attempt-*-sample-audit.json'))
        assert Path(summary['source_audit']).resolve()==ap.resolve() and summary['source_audit_sha256']==sha(ap)
        for p,h in summary['artifacts'].items():assert sha(p)==h
        arrays_path=rp.parent/chain['chain_id']/'states.npz';cp=arrays_path.parent/'coordinates.json'
        assert str(arrays_path) in summary['artifacts'] and str(cp) in summary['artifacts']
        coordinates=json.loads(cp.read_text())
        assert coordinates['input_alignment_sha256']==chain['alignment_sha256']==sha(chain['alignment'])
        assert Path(coordinates['input_alignment']).resolve()==Path(chain['alignment']).resolve()
        with open(chain['alignment']) as handle:
            records=list(SeqIO.parse(handle,'fasta'))
        assert len(records)==len({r.id for r in records})
        assert coordinates['tips']==[dict(tip=r.id,length=len(str(r.seq).replace('-',''))) for r in sorted(records,key=lambda r:r.id)]
        audit=json.loads(ap.read_text());nodes={r['source_node'] for r in audit['candidate_samples']}
        assert coordinates['nodes']==sorted(nodes) and len(nodes)==4
        if common is None:common=coordinates
        assert coordinates==common,'Biological coordinate metadata differs across chains'
        with np.load(arrays_path,allow_pickle=False) as saved:arrays.append(validate_arrays(saved,coordinates,expected))
        assert summary['samples']==len(expected) and summary['state_observations']==arrays[-1].size
        evidence.update({str(p):sha(p) for p in [dp,rp,attempt,arrays_path,cp]})
    return dict(values=np.stack(arrays),coordinates=common,iterations=np.asarray(expected),
        chains=scalar['chains'],evidence=evidence,
        scope='Verified four-chain state input only; no convergence or posterior qualification.')
