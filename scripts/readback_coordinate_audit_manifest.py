#!/usr/bin/env python3
"""Verify the complete coordinate-audit model grid, encoding bytes and summary counts."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['coordinates','mapping','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();rp=a.coordinates/'receipt.json';rh=sha(rp);r=json.loads(rp.read_text());mr=a.mapping/'receipt.json';m=json.loads(mr.read_text())
    if r['status']!='complete_native_3di_coordinate_audit' or r['mapping_receipt_sha256']!=sha(mr):raise ValueError('Incomplete or wrong source mapping')
    provenance=a.mapping/'model_provenance.json'
    if sha(provenance)!=m['artifacts'][provenance.name]:raise ValueError('Mapping provenance changed')
    models=json.loads(provenance.read_text());expected={(x['model_id'],str(x['version'])):x for x in models}
    if len(expected)!=len(models):raise ValueError('Duplicate mapped model')
    summary=a.coordinates/'model_summary.tsv'
    if sha(summary)!=r['artifacts'][summary.name]:raise ValueError('Summary changed')
    seen=set();totals=Counter()
    with summary.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['model_id'],row['version']
            if key in seen or key not in expected:raise ValueError('Wrong/repeated model')
            model=expected[key]
            if row['sequence_sha256']!=model['sequence_sha256'] or int(row['length'])!=model['length']:raise ValueError('Wrong sequence or length')
            if sha(row['encoding_path'])!=row['encoding_sha256']:raise ValueError('Changed encoding bytes')
            values={k:int(row[k]) for k in r['totals']}
            if values['valid_states']+values['invalid_states']!=values['length'] or not 0<=values['valid_feature_plddt70']<=values['valid_focal_plddt70']<=values['valid_states']:raise ValueError('Inconsistent counts')
            totals.update(values);seen.add(key)
    if seen!=set(expected) or len(seen)!=r['models'] or len(seen)!=m['distinct_models'] or dict(totals)!=r['totals'] or sha(rp)!=rh:raise ValueError('Scope/counts or receipt changed')
    result=dict(status='passed_complete_coordinate_audit_model_manifest',models=len(seen),totals=dict(totals),coordinate_receipt_sha256=rh,mapping_receipt_sha256=sha(mr),script_sha256=sha(__file__),scope='Complete model identity/sequence/length grid, every encoding-file hash and all per-model summary sums checked. Does not independently recompute coordinates, native neural states, numerical array contents or PAE qualification.')
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
