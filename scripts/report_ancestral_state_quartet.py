#!/usr/bin/env python3
"""Report all categorical coordinates at both burn-in cutoffs."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import time
import arviz as az
import numpy as np
from ancestral_chain_diagnostics import diagnose,sha
from ancestral_categorical_diagnostics import diagnose_states
from ancestral_state_patterns import group_patterns,reconstruct_patterns


def write_json(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();m=json.loads(args.manifest.read_text());chains=m['chains']
    assert m['status']=='provenance_checked_state_quartet'
    assert len(chains)==len({c['chain_id'] for c in chains})==len({c['seed'] for c in chains})==4
    assert len({c['model_input_identity'] for c in chains})==1
    pins={str(args.manifest):sha(args.manifest),**m['evidence']}
    for p in [__file__,'scripts/ancestral_categorical_diagnostics.py','scripts/ancestral_chain_diagnostics.py','scripts/ancestral_state_patterns.py']:pins[str(p)]=sha(p)
    for p,h in pins.items():assert sha(p)==h
    assert sha(m['arrays'])==m['arrays_sha256'];pins[m['arrays']]=m['arrays_sha256']
    with np.load(m['arrays'],allow_pickle=False) as saved:
        values=saved['values'];iterations=saved['iterations'];free=saved['unanchored_residue_counts']
    coords=m['coordinates'];alphabet=coords['alphabet'];expected=m['expected_iterations']
    assert len(alphabet)==len(set(alphabet)) and np.array_equal(iterations,expected)
    assert expected==sorted(set(expected)) and expected[0]==0
    assert len(coords['nodes'])==len(set(coords['nodes']))==4
    assert values.dtype==np.uint8 and values.shape==(4,len(expected),4,sum(t['length'] for t in coords['tips']))
    assert np.all(values<len(alphabet)) and free.shape==(4,len(expected),4) and np.all(free>=0)
    args.output.mkdir(parents=True,exist_ok=False);outputs={};started=time.monotonic()
    for cutoff in [250,500]:
        assert cutoff<expected[-1]
        retained=values[:,iterations>cutoff];patterns,mapping=group_patterns(retained)
        assert np.array_equal(reconstruct_patterns(patterns,mapping),retained)
        folder=args.output/f'discard-{cutoff}';folder.mkdir()
        np.savez_compressed(folder/'patterns.npz',patterns=patterns,coordinate_pattern_ids=mapping)
        multiplicity=np.bincount(mapping.ravel(),minlength=len(patterns));pattern_counts=Counter();coordinate_counts=Counter()
        with gzip.open(folder/'diagnostics.jsonl.gz','wt') as handle:
            for i,pattern in enumerate(patterns):
                result=diagnose_states(pattern,alphabet);pattern_counts[result['status']]+=1;coordinate_counts[result['status']]+=int(multiplicity[i])
                handle.write(json.dumps(dict(pattern_id=i,coordinate_multiplicity=int(multiplicity[i]),diagnostic=result),allow_nan=False,separators=(',',':'))+'\n')
        # Full output readback: coordinate mapping, report ids, state frequencies
        # and retained denominators; numerical diagnostic fixtures are separate.
        with np.load(folder/'patterns.npz',allow_pickle=False) as saved:
            assert np.array_equal(reconstruct_patterns(saved['patterns'],saved['coordinate_pattern_ids']),retained)
        read_coordinates=0
        with gzip.open(folder/'diagnostics.jsonl.gz','rt') as handle:
            count=0
            for count,line in enumerate(handle,1):
                row=json.loads(line);i=count-1;assert row['pattern_id']==i and row['coordinate_multiplicity']==int(multiplicity[i])
                d=row['diagnostic'];assert d['draws_per_chain']==retained.shape[1]
                for state,label in enumerate(alphabet):
                    counts=np.count_nonzero(patterns[i]==state,axis=1)
                    assert d['indicators'][label]['counts_per_chain']==counts.tolist()
                    assert d['indicators'][label]['frequency_per_chain']==(counts/retained.shape[1]).tolist()
                read_coordinates+=row['coordinate_multiplicity']
            assert count==len(patterns) and read_coordinates==mapping.size
        unanchored={node:diagnose(free[:,iterations>cutoff,n]) for n,node in enumerate(coords['nodes'])}
        write_json(folder/'summary.json',dict(discard_through=cutoff,retained_samples_per_chain=int(retained.shape[1]),patterns=len(patterns),coordinates=int(mapping.size),pattern_status_counts=dict(pattern_counts),coordinate_status_counts=dict(coordinate_counts),unanchored_count_screens=unanchored,scope='Coordinate counts retain duplicate/correlated anchors; not independent biological sites.'))
        outputs[str(cutoff)]=dict(patterns=len(patterns),coordinates=int(mapping.size),retained_samples_per_chain=int(retained.shape[1]))
    for p,h in pins.items():assert sha(p)==h
    write_json(args.output/'receipt.json',dict(status='both_cutoff_categorical_reports_complete_not_posterior_qualification',pins=pins,outputs=outputs,arviz_version=az.__version__,numpy_version=np.__version__,elapsed_seconds=time.monotonic()-started,artifacts={str(p.relative_to(args.output)):sha(p) for p in args.output.rglob('*') if p.is_file()},scope='All declared node/anchor coordinates and states retained. Source-manifest provenance must be independently verified. Marginal diagnostics do not qualify the joint posterior.'))


if __name__=='__main__':main()
