#!/usr/bin/env python3
"""Prepare a provenance-checked quartet and invoke the locked diagnostic environment."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
from ancestral_chain_attempt import sha,write_json
from read_ancestral_state_quartet import read_quartet


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['producer-plan','extraction-plan','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--group',required=True);p.add_argument('--python',required=True);args=p.parse_args()
    producer=json.loads(args.producer_plan.read_text());extraction=json.loads(args.extraction_plan.read_text())
    for plan in [producer,extraction]:
        for name,h in plan['pins'].items():assert sha(name)==h
    assert sha(extraction['producer_plan'])==sha(args.producer_plan)
    state_root=Path(extraction['output']);assert (state_root/'run_plan.json').read_bytes()==args.extraction_plan.read_bytes()
    jobs=[j for j in json.loads(Path(producer['jobs']).read_text()) if j['config']['model_input_identity']==args.group]
    quartet=read_quartet(jobs,producer['output'],state_root,producer['iterations'],sha(producer['mapping']))
    assert quartet is not None,'Four verified extracted chains required'
    free=[]
    for job in jobs:
        chain=job['chain']['chain_id'];d=json.loads((state_root/chain/'disposition.json').read_text())
        path=Path(d['receipt']).parent/chain/'states.npz';assert sha(path)==quartet['evidence'][str(path)]
        with np.load(path,allow_pickle=False) as saved:free.append(saved['unanchored_residue_counts'])
    args.output.mkdir(parents=True,exist_ok=False);arrays=args.output/'quartet.npz'
    np.savez_compressed(arrays,values=quartet['values'],iterations=quartet['iterations'],unanchored_residue_counts=np.stack(free))
    evidence=dict(quartet['evidence'])
    for path in [args.producer_plan,args.extraction_plan,Path(__file__),Path('scripts/read_ancestral_state_quartet.py')]:evidence[str(path)]=sha(path)
    manifest=args.output/'manifest.json';write_json(manifest,dict(status='provenance_checked_state_quartet',chains=quartet['chains'],coordinates=quartet['coordinates'],expected_iterations=quartet['iterations'].tolist(),arrays=str(arrays),arrays_sha256=sha(arrays),evidence=evidence))
    subprocess.run([args.python,'scripts/report_ancestral_state_quartet.py','--manifest',str(manifest),'--output',str(args.output/'reports')],check=True)
    rp=args.output/'reports/receipt.json';report=json.loads(rp.read_text());assert report['status']=='both_cutoff_categorical_reports_complete_not_posterior_qualification'
    for name,h in report['artifacts'].items():assert sha(rp.parent/name)==h
    write_json(args.output/'receipt.json',dict(status='verified_quartet_categorical_reports_complete_not_posterior_qualification',group=args.group,report=str(rp),report_sha256=sha(rp),manifest=str(manifest),manifest_sha256=sha(manifest),scope='Both cutoffs and all coordinate mappings preserved; no global posterior qualification.'))


if __name__=='__main__':main()
