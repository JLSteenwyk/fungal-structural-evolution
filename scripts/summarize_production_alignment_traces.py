#!/usr/bin/env python3
"""Exact extant geometry traces for verified terminal production chains."""
import argparse
import csv
from io import StringIO
import json
from pathlib import Path
import time
from Bio import SeqIO
from ancestral_chain_attempt import sha,write_json
from ancestral_extant_alignment_geometry import project,signature,signature_distance
from readback_independent_baliphy_chain import blocks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    snapshot=json.loads(args.snapshot.read_text());assert snapshot['status']=='verified_frozen_terminal_chain_subset'
    inputs={r['chain_id']:r for r in json.loads(args.inputs.read_text())}
    pins={str(p):sha(p) for p in [args.snapshot,args.inputs,Path(__file__),Path('scripts/ancestral_extant_alignment_geometry.py')]}
    summaries=[];started=time.monotonic()
    for row in snapshot['rows']:
        chain=inputs[row['chain']];assert chain['seed']==row['seed']
        auditpath=Path(row['audit']);assert sha(auditpath)==row['audit_sha256']
        audit=json.loads(auditpath.read_text());rp=Path(audit['attempt_receipt']);assert sha(rp)==audit['attempt_receipt_sha256']
        receipt=json.loads(rp.read_text());assert receipt['exit_code']==0
        for name,digest in receipt['artifacts'].items():assert sha(rp.parent/name)==digest
        assert sha(chain['alignment'])==chain['alignment_sha256']
        records=list(SeqIO.parse(chain['alignment'],'fasta'));original={r.id:str(r.seq).upper() for r in records};assert len(original)==len(records)
        observed={t:s.replace('-','') for t,s in original.items()}
        ref=signature(project(original,observed));previous=None;data=[]
        paths=list(rp.parent.glob('independent-chain-*/C1.P1.fastas'));assert len(paths)==1
        for iteration,text in blocks(paths[0]):
            records=list(SeqIO.parse(StringIO(text.lstrip()),'fasta'));sequences={r.id:str(r.seq).upper() for r in records};assert len(sequences)==len(records)
            geometry=project(sequences,observed);current=signature(geometry)
            data.append(dict(iteration=iteration,distance_from_input=signature_distance(ref,current),
                distance_from_previous_saved='' if previous is None else signature_distance(previous,current),
                extant_alignment_columns=len(next(iter(geometry.values())))))
            previous=current
        assert [r['iteration'] for r in data]==list(range(0,audit['iterations']+1,10))
        assert len(data)==row['saved_alignments']
        p=args.output/(row['chain']+'.tsv')
        with p.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows(data)
        summaries.append(dict(chain=row['chain'],family=chain['family'],input_group=chain['effective_input_group'],
            original_configuration_ids=chain['original_configuration_ids'],seed=chain['seed'],tips=len(observed),
            residues=sum(map(len,observed.values())),saved_samples=len(data),
            initial_input_distance=data[0]['distance_from_input'],final_input_distance=data[-1]['distance_from_input'],
            minimum_input_distance=min(r['distance_from_input'] for r in data),maximum_input_distance=max(r['distance_from_input'] for r in data),
            unchanged_successive_saved_pairs=sum(r['distance_from_previous_saved']==0 for r in data[1:]),
            compared_successive_saved_pairs=len(data)-1,trace=str(p),trace_sha256=sha(p),
            source_audit=str(auditpath),source_audit_sha256=sha(auditpath),source_samples=str(paths[0]),source_samples_sha256=sha(paths[0])))
        print(json.dumps({'chains':len(summaries),'elapsed_seconds':time.monotonic()-started}),flush=True)
    for p,h in pins.items():assert sha(p)==h
    write_json(args.output/'receipt.json',dict(status='production_geometry_traces_complete_not_convergence',
        pins=pins,chains=len(summaries),saved_samples=sum(r['saved_samples'] for r in summaries),summaries=summaries,
        elapsed_seconds=time.monotonic()-started,
        scope='Exact integer geometry distances to the fixed observed input alignment and previous saved alignment; '
        'all samples include burn-in. Nonzero change demonstrates movement only. No independent-chain agreement, '
        'stationarity, effective sample size, or ancestral-state mixing conclusion.'))


if __name__=='__main__':main()
