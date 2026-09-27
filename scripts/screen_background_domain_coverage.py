#!/usr/bin/env python3
"""Preserve both orders and screen domain coverage; never certify diagnostic data for inference."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def reasons(statuses, lengths, original_lengths, retained_lengths, minimum, fraction):
    """Both orders must satisfy length and coverage relative to ORIGINAL intervals."""
    result=[]
    if any(s=='input_unavailable' for s in statuses):result.append('input_unavailable')
    if any(s=='outside_printed_rounding' for s in statuses):result.append('numeric_discrepancy')
    if any(s not in {'within_printed_rounding','outside_printed_rounding','input_unavailable'} for s in statuses):
        raise ValueError('Unknown disposition')
    if len(lengths)==2:
        if min(lengths)<minimum:result.append('short_alignment')
        cutoff=Fraction(str(fraction))
        if any(n*cutoff.denominator < cutoff.numerator*d for n in lengths for d in original_lengths):result.append('low_original_interval_coverage')
        if any(n>min(retained_lengths) for n in lengths):raise ValueError('Alignment exceeds retained input')
    elif not result:raise ValueError('Missing order')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed input: '+path)
    verify()
    diagnostic=Path(plan['diagnostic']);dr=json.loads((diagnostic/'receipt.json').read_text())
    if dr['status']!='complete_background_domain_rmsd_diagnostic_not_scientific_acceptance' or dr['scientific_eligibility'] is not False:raise ValueError('Wrong diagnostic source')
    if sha(diagnostic/'numeric_readback.tsv')!=dr['artifacts']['numeric_readback.tsv']:raise ValueError('Changed numeric table')
    dp=json.loads(Path(plan['diagnostic_plan']).read_text());sp=json.loads(Path(dp['source_plan']).read_text())
    if dr['plan_sha256']!=sha(plan['diagnostic_plan']) or dp['output']!=str(diagnostic) or sp['inventory']!=plan['inventory'] or sp['inputs']!=plan['inputs']:raise ValueError('Diagnostic input provenance differs')
    inventory=Path(plan['inventory']);ir=json.loads((inventory/'receipt.json').read_text())
    for name in ['intervals.jsonl','domain_pairs.tsv']:
        if sha(inventory/name)!=ir['artifacts'][name]:raise ValueError('Changed inventory')
    intervals={}
    for line in (inventory/'intervals.jsonl').open():
        r=json.loads(line);key=r['interval_id']
        if key in intervals or r['length']!=r['end']-r['start']+1:raise ValueError('Invalid interval')
        intervals[key]=r
    material=Path(plan['inputs']);mr=json.loads((material/'receipt.json').read_text())
    if sha(material/'inputs.jsonl')!=mr['artifacts']['inputs.jsonl']:raise ValueError('Changed inputs')
    inputs={}
    for line in (material/'inputs.jsonl').open():
        r=json.loads(line);key=r['interval_id'],r['mask']
        if key in inputs or r['interval_length']!=intervals[key[0]]['length']:raise ValueError('Invalid input')
        inputs[key]={'status':r['status'],'retained_residues':r['retained_residues']}
    nums={}
    with (diagnostic/'numeric_readback.tsv').open() as f:
        reader=csv.DictReader(f,delimiter='\t');metrics=[c for c in reader.fieldnames if c not in ['pair_key','mask','order']]
        for r in reader:
            key=r['pair_key'],r['mask'],int(r['order'])
            if key in nums:raise ValueError('Repeated numeric record')
            nums[key]=r
    pairs=list(csv.DictReader((inventory/'domain_pairs.tsv').open(),delimiter='\t'))
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    fields=['pair_key','mask','interval_a','interval_b','length_a','length_b','retained_a','retained_b','input_status_a','input_status_b']
    for order in [0,1]:fields += [f'order{order}_status']+[f'order{order}_{m}' for m in metrics if m!='rmsd_status']+[f'order{order}_original_coverage_a',f'order{order}_original_coverage_b']
    for spec in plan['screens']:fields += [spec['id']+'_pass',spec['id']+'_exclusions']
    seen=set();counts=Counter();exclusions=Counter();pair_flags={};rows=0;unavailable=0;quarantined=0
    with (out/'pair_mask_coverage.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for pair in pairs:
            key=pair['domain_pair_key'];ends=[pair['interval_a'],pair['interval_b']]
            if key!=hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest():raise ValueError('Invalid pair identity')
            if key in pair_flags:raise ValueError('Repeated pair')
            pair_flags[key]={}
            for mask in ['full','plddt70']:
                lens=[intervals[e]['length'] for e in ends];inps=[inputs[e,mask] for e in ends];retained=[x['retained_residues'] for x in inps]
                ready=all(x['status']=='ready' for x in inps)
                row=dict(pair_key=key,mask=mask,interval_a=ends[0],interval_b=ends[1],length_a=lens[0],length_b=lens[1],retained_a=retained[0],retained_b=retained[1],input_status_a=inps[0]['status'],input_status_b=inps[1]['status'])
                statuses=[];lengths=[]
                for order in [0,1]:
                    nk=key,mask,order;n=nums.get(nk)
                    if bool(n)!=ready:raise ValueError('Missing or unexpected numeric disposition')
                    status=n['rmsd_status'] if n else 'input_unavailable';statuses.append(status);row[f'order{order}_status']=status
                    if n:
                        seen.add(nk);length=int(n['aligned_length']);lengths.append(length)
                        for m in metrics:
                            if m!='rmsd_status':row[f'order{order}_{m}']=n[m]
                        row[f'order{order}_original_coverage_a']=length/lens[0];row[f'order{order}_original_coverage_b']=length/lens[1]
                unavailable+=not ready;quarantined+=('outside_printed_rounding' in statuses)
                pair_flags[key][mask]={}
                for spec in plan['screens']:
                    why=reasons(statuses,lengths,lens,retained,spec['minimum_aligned_residues'],spec['minimum_original_coverage'])
                    passed=not why;row[spec['id']+'_pass']=int(passed);row[spec['id']+'_exclusions']=';'.join(why)
                    counts[mask+':'+spec['id']]+=passed
                    for reason in why:exclusions[mask+':'+spec['id']+':'+reason]+=1
                    pair_flags[key][mask][spec['id']]=passed
                writer.writerow(row);rows+=1
    if seen!=set(nums) or len(seen)!=dr['numerically_checked_alignments'] or rows*2!=dr['directed_dispositions']:raise ValueError('Incomplete cohort')
    joint={s['id']:sum(v['full'][s['id']] and v['plddt70'][s['id']] for v in pair_flags.values()) for s in plan['screens']}
    verify()
    receipt=dict(status='complete_background_domain_coverage_screen_pending_independent_readback',plan_sha256=ph,diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),pair_mask_rows=rows,distinct_pairs=len(pairs),unavailable_pair_masks=unavailable,quarantined_pair_masks=quarantined,pass_counts=dict(counts),exclusion_counts=dict(exclusions),both_masks_pass_counts=joint,artifacts={'pair_mask_coverage.tsv':sha(out/'pair_mask_coverage.tsv')},biological_inference_eligible=False,scope='Descriptive coverage screen only. Original diagnostic remains non-passing; mismatch quarantined. Both orders required. Coverage uses original interval length, not confidence-retained length. No order selection, native metric correction, domain boundary validation, geometric-rank test or duplication effect inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
