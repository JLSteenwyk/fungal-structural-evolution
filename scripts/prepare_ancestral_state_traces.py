#!/usr/bin/env python3
"""Preserve residue-anchored categorical traces; no convergence qualification."""
import argparse
from io import StringIO
import json
from pathlib import Path
import time
import numpy as np
from Bio import SeqIO
from ancestral_chain_attempt import sha,write_json
from ancestral_residue_anchors import anchored_columns
from readback_independent_baliphy_chain import blocks

ALPHABET='ACDEFGHIKLMNPQRSTVWYX-'


def project_states(sequences, observed, candidates):
    tips=sorted(observed);nodes=sorted(candidates)
    offsets={};size=0
    for tip in tips:offsets[tip]=size;size+=len(observed[tip])
    columns=anchored_columns(sequences,observed,candidates)
    result=np.full((len(nodes),size),255,dtype=np.uint8)
    unanchored=np.zeros(len(nodes),dtype=np.int32)
    for column in columns:
        if column['unanchored']:
            for n,node in enumerate(nodes):unanchored[n]+=column['states'][node]!='-'
        else:
            indices=[offsets[t]+position-1 for t,position in column['anchors']]
            for n,node in enumerate(nodes):result[n,indices]=ALPHABET.index(column['states'][node])
    assert not np.any(result==255)
    # Independent traversal selects ancestral letters at each tip's nongap
    # columns. It neither calls the anchor builder nor reuses its positions.
    lookup=np.full(256,255,dtype=np.uint8)
    for i,a in enumerate(ALPHABET):lookup[ord(a)]=i
    for n,node in enumerate(nodes):
        sampled=np.frombuffer(sequences[candidates[node]].encode('ascii'),dtype=np.uint8)
        direct=np.concatenate([lookup[sampled[np.frombuffer(sequences[t].encode('ascii'),dtype=np.uint8)!=ord('-')]] for t in tips])
        assert np.array_equal(result[n],direct)
        reconstructed=''.join(c['states'][node] for c in columns)
        assert reconstructed==sequences[candidates[node]]
        extant_present=np.any(np.array([list(sequences[t]) for t in tips])!='-',axis=0)
        assert unanchored[n]==int(np.count_nonzero((sampled!=ord('-')) & ~extant_present))
    return result,unanchored


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    snapshot=json.loads(args.snapshot.read_text());assert snapshot['status']=='verified_frozen_terminal_chain_subset'
    inputs={r['chain_id']:r for r in json.loads(args.inputs.read_text())}
    pins={str(p):sha(p) for p in [args.snapshot,args.inputs,Path(__file__),Path('scripts/ancestral_residue_anchors.py')]}
    summaries=[];started=time.monotonic()
    for row in snapshot['rows']:
        chain=inputs[row['chain']];assert chain['seed']==row['seed']
        ap=Path(row['audit']);assert sha(ap)==row['audit_sha256'];audit=json.loads(ap.read_text())
        rp=Path(audit['attempt_receipt']);assert sha(rp)==audit['attempt_receipt_sha256']
        receipt=json.loads(rp.read_text());assert receipt['exit_code']==0
        for name,h in receipt['artifacts'].items():assert sha(rp.parent/name)==h
        assert sha(chain['alignment'])==chain['alignment_sha256']
        observed={r.id:str(r.seq).replace('-','').upper() for r in SeqIO.parse(chain['alignment'],'fasta')}
        by_iteration={}
        for r in audit['candidate_samples']:by_iteration.setdefault(r['iteration'],{})[r['source_node']]=r['runtime_node']
        paths=list(rp.parent.glob('independent-chain-*/C1.P1.fastas'));assert len(paths)==1
        traces=[];unanchored=[];iterations=[];candidates=None
        for iteration,text in blocks(paths[0]):
            records=list(SeqIO.parse(StringIO(text.lstrip()),'fasta'));sequences={r.id:str(r.seq).upper() for r in records};assert len(records)==len(sequences)
            current=by_iteration[iteration];assert len(current)==4
            if candidates is None:candidates=current
            assert current==candidates
            values,free=project_states(sequences,observed,candidates)
            traces.append(values);unanchored.append(free);iterations.append(iteration)
        assert iterations==list(range(0,audit['iterations']+1,10))
        values=np.stack(traces);free=np.stack(unanchored);its=np.asarray(iterations)
        counts={};variable={}
        for cutoff in [250,500]:
            retained=values[its>cutoff]
            c=np.stack([np.count_nonzero(retained==state,axis=0) for state in range(len(ALPHABET))],axis=-1).astype(np.uint16)
            assert np.all(c.sum(axis=-1)==len(retained))
            counts[f'counts_after_{cutoff}']=c
            variable[str(cutoff)]=int(np.count_nonzero(np.count_nonzero(c,axis=-1)>1))
        target=args.output/row['chain'];target.mkdir()
        np.savez_compressed(target/'states.npz',states=values,iterations=its,unanchored_residue_counts=free,**counts)
        # Serialization roundtrip includes the complete temporal trace.
        with np.load(target/'states.npz',allow_pickle=False) as saved:
            assert np.array_equal(saved['states'],values) and np.array_equal(saved['iterations'],its)
            assert np.array_equal(saved['unanchored_residue_counts'],free)
            for name,c in counts.items():assert np.array_equal(saved[name],c)
        coordinate=dict(alphabet=ALPHABET,nodes=sorted(candidates),tips=[dict(tip=t,length=len(observed[t])) for t in sorted(observed)],
            ordering='state axes: saved iteration, source node, concatenated ungapped input residues in listed tip order; positions are one-based within each tip',
            input_alignment=chain['alignment'],input_alignment_sha256=chain['alignment_sha256'])
        write_json(target/'coordinates.json',coordinate)
        summaries.append(dict(chain=row['chain'],samples=len(iterations),candidate_anchor_coordinates=int(values.shape[1]*values.shape[2]),
            state_observations=int(values.size),variable_coordinates=variable,unanchored_residue_observations=int(free.sum()),
            source_audit=str(ap),source_audit_sha256=sha(ap),
            artifacts={str(p):sha(p) for p in [target/'states.npz',target/'coordinates.json']}))
        print(json.dumps({'chains':len(summaries),'elapsed_seconds':time.monotonic()-started}),flush=True)
    for p,h in pins.items():assert sha(p)==h
    write_json(args.output/'receipt.json',dict(status='complete_independently_checked_anchored_state_traces',pins=pins,
        chains=len(summaries),state_observations=sum(r['state_observations'] for r in summaries),summaries=summaries,
        elapsed_seconds=time.monotonic()-started,
        scope='Trace of candidate state at each extant residue anchor, including gap and X, with both burn-in count summaries. '
        'Anchors can duplicate the same ancestral residue and are not independent sites. Unanchored residues have counts only, '
        'no cross-sample identity; full original alignments remain authoritative. No convergence or posterior qualification.'))


if __name__=='__main__':main()
