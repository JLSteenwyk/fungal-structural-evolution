#!/usr/bin/env python3
"""Verify every domain boundary link, interval descriptor and deduplicated interval pair."""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed readback plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed pinned source: '+path)
    verify();sp=json.loads(Path(plan['source_plan']).read_text());folder=Path(sp['output']);r=json.loads((folder/'receipt.json').read_text());controls=Path(sp['controls']);cp=json.loads(Path(sp['control_plan']).read_text())
    if r['status']!='complete_duplication_domain_pair_inventory_pending_readback' or r['plan_sha256']!=sha(plan['source_plan']) or r['source_control_receipt_sha256']!=sha(controls/'receipt.json'):raise ValueError('Unbound inventory')
    audit=json.loads(Path(sp['readback']).read_text())
    if audit['status']!='passed_full_duplication_domain_control_readback' or audit['producer_receipt_sha256']!=sha(controls/'receipt.json'):raise ValueError('Domain controls not audited')
    for name,h in r['artifacts'].items():
        if sha(folder/name)!=h:raise ValueError('Changed inventory artifact')
    models={}
    for source in [cp['base_queue'],cp['inventory']]:
        with (Path(source)/'models.jsonl').open() as handle:
            for line in handle:
                m=json.loads(line);key=m['model_id'],m['version']
                if key in models and models[key]!=m:raise ValueError('Shared model differs')
                models[key]=m
    intervals={}
    with (folder/'intervals.jsonl').open() as handle:
        for line in handle:
            row=json.loads(line);key=row['model_id'],row['version'];m=models[key];start,end=row['start'],row['end']
            if not 1<=start<=end<=m['length']:raise ValueError('Invalid interval bounds')
            ident=hashlib.sha256(json.dumps([key[0],key[1],m['sequence_sha256'],start,end],separators=(',',':')).encode()).hexdigest()
            expected=dict(interval_id=ident,model_id=key[0],version=key[1],sequence_sha256=m['sequence_sha256'],source_sha256=m['sha256'],source_path=m['path'],original_length=m['length'],start=start,end=end,length=end-start+1)
            if row!=expected or ident in intervals:raise ValueError('Interval identity/provenance differs')
            intervals[ident]=row
    annotations={}
    with (controls/'model_annotations.jsonl').open() as handle:
        for line in handle:
            row=json.loads(line);key=row['model_id'],row['version']
            annotations[key]={policy:{h['hit_id']:h for h in hits} for policy,hits in row['policies'].items()}
    used=set();expected_pairs={};counts=Counter();source_rows=0
    with (controls/'single_copy_domain_pairs.tsv').open() as source,(folder/'domain_pair_links.tsv').open() as links:
        exported=iter(csv.DictReader(links,delimiter='\t'))
        for row in csv.DictReader(source,delimiter='\t'):
            source_rows+=1
            for boundary in ['alignment','envelope']:
                link=next(exported,None)
                if link is None or any(link[k]!=v for k,v in row.items()) or link['boundary']!=boundary:raise ValueError('Original domain-match link differs')
                ids=[]
                for side in ['left','right']:
                    key=row['model_'+side],int(row['version_'+side]);hit=annotations[key][row['policy']][row['hit_'+side]]
                    if hit['pfam_accession']!=row['pfam_accession'] or not hit['domain_interval_candidate']:raise ValueError('Source hit eligibility differs')
                    ident=link['interval_'+side];interval=intervals[ident]
                    if (interval['model_id'],interval['version'],interval['start'],interval['end'])!=(key[0],key[1],hit[boundary+'_start'],hit[boundary+'_end']):raise ValueError('Boundary interval differs from source hit')
                    used.add(ident);ids.append(ident)
                if ids[0]==ids[1]:raise ValueError('Identical interval pair')
                endpoints=sorted(ids);pair=hashlib.sha256(json.dumps(endpoints,separators=(',',':')).encode()).hexdigest()
                if link['domain_pair_key']!=pair:raise ValueError('Domain pair hash differs')
                if pair in expected_pairs and expected_pairs[pair]!=endpoints:raise ValueError('Pair hash collision')
                expected_pairs[pair]=endpoints;counts[row['policy']+':'+boundary]+=1
        if next(exported,None) is not None:raise ValueError('Extra domain link')
    actual_pairs={}
    with (folder/'domain_pairs.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=row['domain_pair_key']
            if key in actual_pairs:raise ValueError('Duplicate interval pair')
            actual_pairs[key]=[row['interval_a'],row['interval_b']]
    if actual_pairs!=expected_pairs or used!=set(intervals):raise ValueError('Incomplete deduplicated interval/pair universe')
    totals=dict(source_domain_match_rows=source_rows,policy_boundary_links=sum(counts.values()),link_counts=dict(counts),unique_intervals=len(intervals),unique_interval_pairs=len(actual_pairs),unique_models=len({(v['model_id'],v['version']) for v in intervals.values()}),interval_residues=sum(v['length'] for v in intervals.values()),maximum_interval_length=max(v['length'] for v in intervals.values()),planned_directed_mask_dispositions=4*len(actual_pairs))
    if any(r[k]!=v for k,v in totals.items()):raise ValueError('Inventory summary differs')
    verify()
    result=dict(status='passed_full_duplication_domain_pair_inventory_readback',plan_sha256=ph,producer_receipt_sha256=sha(folder/'receipt.json'),**totals,scope='All source Domain-match fields and both boundary expansions reconstructed; every interval bounds/identity/provenance field and exact deduplicated interval/pair sets checked. No domain coordinate extraction, structural alignment, PAE qualification or biological inference.')
    with Path(plan['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in totals.items() if k!='link_counts'}),flush=True)


if __name__=='__main__':main()
