#!/usr/bin/env python3
"""Independently check every domain coverage screen row against numeric and interval sources."""
import argparse,csv,hashlib,json
from collections import Counter
from fractions import Fraction
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_domain_coverage_screen_pending_independent_readback'
    assert receipt['plan_sha256']==sha(args.plan) and receipt['biological_inference_eligible'] is False
    for path,h in plan['pins'].items():assert sha(path)==h,path
    for name,h in receipt['artifacts'].items():assert sha(root/name)==h
    dr=Path(plan['diagnostic'])/'receipt.json';assert sha(dr)==receipt['diagnostic_receipt_sha256']
    pairs={r['domain_pair_key']:(r['interval_a'],r['interval_b']) for r in csv.DictReader((Path(plan['inventory'])/'domain_pairs.tsv').open(),delimiter='\t')}
    lengths={}
    for line in (Path(plan['inventory'])/'intervals.jsonl').open():
        r=json.loads(line);lengths[r['interval_id']]=r['end']-r['start']+1
    inputs={}
    for line in (Path(plan['inputs'])/'inputs.jsonl').open():
        r=json.loads(line);inputs[r['interval_id'],r['mask']]=(r['status'],len(r['original_positions']))
    numbers={}
    for r in csv.DictReader((Path(plan['diagnostic'])/'numeric_readback.tsv').open(),delimiter='\t'):
        numbers[r['pair_key'],r['mask'],int(r['order'])]=r
    seen=set();numeric_seen=set();counts=Counter();exclusions=Counter();passed_masks={s['id']:Counter() for s in plan['screens']};unavailable=quarantine=0
    for row in csv.DictReader((root/'pair_mask_coverage.tsv').open(),delimiter='\t'):
        pk=row['pair_key'];mask=row['mask'];key=pk,mask
        assert mask in ['full','plddt70'] and key not in seen;seen.add(key)
        a,b=pairs[pk];assert [row['interval_a'],row['interval_b']]==[a,b]
        for suffix,e in zip('ab',[a,b]):
            assert int(row['length_'+suffix])==lengths[e]
            assert row['input_status_'+suffix]==inputs[e,mask][0]
            assert int(row['retained_'+suffix])==inputs[e,mask][1]
        ready=inputs[a,mask][0]==inputs[b,mask][0]=='ready'
        records=[numbers.get((pk,mask,o)) for o in [0,1]]
        assert all((x is not None)==ready for x in records)
        bad=False;ns=[]
        for order,record in enumerate(records):
            prefix=f'order{order}_'
            if record is None:
                assert row[prefix+'status']=='input_unavailable'
                assert all(not value for col,value in row.items() if col.startswith(prefix) and col!=prefix+'status')
                continue
            numeric_seen.add((pk,mask,order));assert row[prefix+'status']==record['rmsd_status']
            bad|=record['rmsd_status']=='outside_printed_rounding'
            for col,value in record.items():
                if col not in ['pair_key','mask','order','rmsd_status']:assert row[prefix+col]==value
            n=int(record['aligned_length']);ns.append(n)
            for suffix,e in zip('ab',[a,b]):assert abs(float(row[prefix+'original_coverage_'+suffix])-n/lengths[e])<1e-14
        unavailable+=not ready;quarantine+=bad
        for spec in plan['screens']:
            why=[]
            if not ready:why.append('input_unavailable')
            if bad:why.append('numeric_discrepancy')
            if ready:
                if any(n<spec['minimum_aligned_residues'] for n in ns):why.append('short_alignment')
                cutoff=Fraction(str(spec['minimum_original_coverage']))
                if min(Fraction(n,lengths[e]) for n in ns for e in [a,b])<cutoff:why.append('low_original_interval_coverage')
            sid=spec['id'];passed=not why
            assert row[sid+'_pass']==str(int(passed)) and row[sid+'_exclusions']==';'.join(why)
            counts[mask+':'+sid]+=passed
            for reason in why:exclusions[mask+':'+sid+':'+reason]+=1
            if passed:passed_masks[sid][pk]+=1
    assert seen=={(pk,m) for pk in pairs for m in ['full','plddt70']}
    assert numeric_seen==set(numbers)
    assert len(seen)==receipt['pair_mask_rows'] and len(pairs)==receipt['distinct_pairs']
    assert unavailable==receipt['unavailable_pair_masks'] and quarantine==receipt['quarantined_pair_masks']
    assert dict(counts)==receipt['pass_counts'] and dict(exclusions)==receipt['exclusion_counts']
    both={sid:sum(n==2 for n in counter.values()) for sid,counter in passed_masks.items()}
    assert both==receipt['both_masks_pass_counts']
    result=dict(status='passed_complete_domain_coverage_screen_readback',producer_receipt_sha256=sha(root/'receipt.json'),plan_sha256=sha(args.plan),checker_sha256=sha(__file__),rows_checked=len(seen),numeric_rows_checked=len(numeric_seen),unavailable_pair_masks=unavailable,quarantined_pair_masks=quarantine,pass_counts=dict(counts),both_masks_pass_counts=both,biological_inference_eligible=False,scope='Every pair/mask/order, copied numeric field, original interval denominator, mask count, rational threshold decision, exclusion reason and summary checked independently. Does not recalculate native structures or confer scientific eligibility.')
    with args.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
