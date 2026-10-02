#!/usr/bin/env python3
"""Reconstruct all full sequence sets and every original triad linkage."""
import argparse
import hashlib
import json
from pathlib import Path

from full_triad_sequence_sources import load_original
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def inspect(plan_path, output):
    plan=json.loads(Path(plan_path).read_text())
    catalog,ready,inputs,mapping,bindings=load_original(plan,plan_path)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_triad_sequence_catalog_pending_independent_readback'
    assert receipt['plan_sha256']==sha(plan_path)
    for name,digest in receipt['artifacts'].items(): bind(bindings,root/name,digest)
    expected={}; scheduled=0; lengths=[]; letters=set()
    with (root/'triad_sequence_links.jsonl').open() as f:
        for source in catalog:
            line=next(f,None);assert line is not None
            row=json.loads(line)
            assert all(row[k]==v for k,v in source.items())
            if not source['source_design_ready_links']:
                assert row['sequence_control_disposition']=='no_source_ready_logical_link'
                assert row['sequence_set_id'] is None and row['role_to_sequence_indices'] is None
                continue
            ordered=sorted(tuple(m) for m in source['models'])
            sid=hashlib.sha256(json.dumps(ordered,separators=(',',':')).encode()).hexdigest()
            seq=[inputs[(*m,'full')]['sequence'] for m in ordered]
            expected[sid]=dict(sequence_set_id=sid,models=[list(m) for m in ordered],sequences=seq,
                sequence_sha256=[hashlib.sha256(s.encode()).hexdigest() for s in seq],original_lengths=[len(s) for s in seq])
            assert row['sequence_control_disposition']=='scheduled_full_sequence_alignment'
            assert row['sequence_set_id']==sid
            assert [ordered[i] for i in row['role_to_sequence_indices']]==[tuple(m) for m in source['models']]
            assert sorted(row['role_to_sequence_indices'])==[0,1,2]
            scheduled+=1
        assert next(f,None) is None
    with (root/'sequence_sets.jsonl').open() as f:
        for sid in sorted(expected):
            line=next(f,None);assert line is not None
            row=json.loads(line);assert row==expected[sid]
            lengths.append(row['original_lengths']);letters.update(''.join(row['sequences']))
        assert next(f,None) is None
    summary=dict(ordered_model_triads=len(catalog),source_ready_triads=scheduled,
        unscheduled_original_triads=len(catalog)-scheduled,unique_sequence_model_sets=len(expected),
        native_alignment_states=12*len(expected),target_contexts=mapping['target_contexts'],
        reference_tie_records=mapping['reference_tie_records'],duplicate_reference_links=mapping['duplicate_reference_links'],
        maximum_sequence_length=max(max(v) for v in lengths),maximum_three_sequence_length=max(map(sum,lengths)),
        sequence_letters=''.join(sorted(letters)))
    assert all(receipt[k]==v for k,v in summary.items())
    bind(bindings,rp);verify(bindings)
    result=dict(status='passed_full_triad_sequence_catalog_readback',**summary,
        plan_sha256=sha(plan_path),producer_receipt_sha256=sha(rp),source_hashes=bindings,
        scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();inspect(a.plan,a.output)
