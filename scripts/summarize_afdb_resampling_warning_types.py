#!/usr/bin/env python3
"""Classify all retained native warning lines, preserving unknown warnings."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re
from review_afdb_paired_resampling_outputs import sha


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--review',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads((args.review/'receipt.json').read_text())
    assert receipt['status']=='completed_interval_readback_and_native_warning_census'
    for name,h in receipt['artifacts'].items():assert sha(args.review/name)==h
    patterns={
        'near_zero_internal_branches':r'WARNING: \d+ near-zero internal branches \(<[^)]+\) should be treated with caution',
        'rare_alignment_states':r'WARNING: States\(s\) .+ rarely appear in alignment and may cause numerical problems',
        'gap_or_ambiguity_over_50_percent':r'WARNING: \d+ sequences contain more than 50% gaps/ambiguity',
        'long_branches_over_9_8':r'WARNING: \d+ too long branches \(>9\.8000\) should be treated with caution!'}
    counts=Counter();byfit=Counter();markers={k:set() for k in patterns};unknown=Counter();rows=warned=0
    seen=set()
    with gzip.open(args.review/'fit_warning_inventory.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line);key=(r['marker'],r['block_length'],r['replicate'],r['fit']);assert key not in seen;seen.add(key)
            rows+=1;warned+=bool(r['warnings']);kinds=set()
            for warning in r['warnings']:
                matches=[k for k,p in patterns.items() if re.fullmatch(p,warning)]
                assert len(matches)<=1
                if not matches:unknown[warning]+=1
                else:kinds.add(matches[0])
            for kind in kinds:
                counts[kind]+=1;byfit[r['fit']+'|'+kind]+=1;markers[kind].add(r['marker'])
    assert rows==receipt['fits'] and warned==receipt['fits_with_warnings']
    result=dict(status='complete_native_warning_classification' if not unknown else 'unknown_native_warnings_require_review',
        fits=rows,fits_with_warnings=warned,counts=dict(counts),by_fit=dict(byfit),
        markers={k:sorted(v) for k,v in markers.items()},unknown_warnings=dict(unknown),patterns=patterns,
        source_receipt_sha256=sha(args.review/'receipt.json'),source_inventory_sha256=sha(args.review/'fit_warning_inventory.jsonl.gz'),
        script_sha256=sha(__file__),scope='Counts are fits containing each warning class; classes overlap. Warnings remain limitations, not automatically excluded fits or proof of failed optimization. Fixed-topology conditional intervals require additional evolutionary model and branch-information checks.')
    assert not args.output.exists()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','fits','counts','by_fit','unknown_warnings']},indent=2))


if __name__=='__main__':main()
