#!/usr/bin/env python3
"""Reconstruct exact diagnostic inputs with bounded standard-library streaming."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as f:
        while part:=f.read(1024*1024):result.update(part)
    return result.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    spec=json.loads(args.manifest.read_text());source=Path(spec['original_audit_export'])
    assert sha(source)==spec['original_audit_export_sha256']
    expected=spec['candidate_audits'];lookup={r['audit_id']:r for r in expected}
    assert len(lookup)==len(expected)==10
    recovered=[]
    with gzip.open(source,'rt') as f:
        for line in f:
            record=json.loads(line)
            if record['audit_id'] in lookup:
                encoded=json.dumps(record,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
                assert hashlib.sha256(encoded).hexdigest()==lookup[record['audit_id']]['record_sha256']
                assert record['cohort_id']==spec['cohort_id'] and record['tree']=='mafft_guide'
                recovered.append(record)
                if len(recovered)==len(expected):break
    assert [r['audit_id'] for r in recovered]==[r['audit_id'] for r in expected]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(recovered,f,indent=2);f.write('\n')
    assert sha(args.output)==spec['candidate_json_sha256']
    print(json.dumps(dict(status='reconstructed_exact_original_diagnostic_candidate_bytes',
        candidates=len(recovered),output=str(args.output),sha256=sha(args.output),scientific_eligibility=False)))


if __name__=='__main__':main()
