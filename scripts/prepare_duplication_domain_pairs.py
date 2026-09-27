#!/usr/bin/env python3
"""Expand verified domain matches to alignment/envelope interval comparisons without losing policy links."""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def interval(model,start,end):
    if not 1<=start<=end<=model['length']:raise ValueError('Out-of-range domain interval')
    key=digest([model['model_id'],model['version'],model['sequence_sha256'],start,end])
    return dict(interval_id=key,model_id=model['model_id'],version=model['version'],sequence_sha256=model['sequence_sha256'],source_sha256=model['sha256'],source_path=model['path'],original_length=model['length'],start=start,end=end,length=end-start+1)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        if sha(a.plan)!=ph:raise ValueError('Changed plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();source=Path(plan['controls']);r=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text())
    if audit['status']!='passed_full_duplication_domain_control_readback' or audit['producer_receipt_sha256']!=sha(source/'receipt.json'):raise ValueError('Unbound domain controls')
    for name,h in r['artifacts'].items():
        if sha(source/name)!=h:raise ValueError('Changed domain control artifact')
    cp=json.loads(Path(plan['control_plan']).read_text())
    if r['plan_sha256']!=sha(plan['control_plan']):raise ValueError('Domain control plan differs')
    models={}
    for root in [Path(cp['base_queue']),Path(cp['inventory'])]:
        receipt=json.loads((root/'receipt.json').read_text());path=root/'models.jsonl'
        if sha(path)!=receipt['artifacts']['models.jsonl']:raise ValueError('Changed model descriptors')
        with path.open() as f:
            for line in f:
                m=json.loads(line);key=m['model_id'],m['version']
                if key in models and models[key]!=m:raise ValueError('Inconsistent model identity')
                models[key]=m
    annotations={}
    with (source/'model_annotations.jsonl').open() as f:
        for line in f:
            row=json.loads(line);key=row['model_id'],row['version']
            if key in annotations or models[key]['sequence_sha256']!=row['sequence_sha256']:raise ValueError('Model annotation binding differs')
            annotations[key]={p:{h['hit_id']:h for h in hits} for p,hits in row['policies'].items()}
    if set(models)!=set(annotations):raise ValueError('Incomplete annotation model universe')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);intervals={};pairs={};counts=Counter();source_rows=0
    with (source/'single_copy_domain_pairs.tsv').open() as f,(out/'domain_pair_links.tsv').open('w') as target:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames+['boundary','interval_left','interval_right','domain_pair_key']
        writer=csv.DictWriter(target,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in reader:
            source_rows+=1;keys=[(row['model_'+s],int(row['version_'+s])) for s in ['left','right']]
            hits=[annotations[key][row['policy']][row['hit_'+s]] for key,s in zip(keys,['left','right'])]
            for h in hits:
                if h['pfam_accession']!=row['pfam_accession'] or not h['domain_interval_candidate']:raise ValueError('Domain candidate binding differs')
            for boundary in ['alignment','envelope']:
                ids=[]
                for key,h in zip(keys,hits):
                    record=interval(models[key],h[boundary+'_start'],h[boundary+'_end']);iid=record['interval_id']
                    if iid in intervals and intervals[iid]!=record:raise ValueError('Interval hash collision')
                    intervals[iid]=record;ids.append(iid)
                if ids[0]==ids[1]:raise ValueError('Unexpected identical model interval pair')
                endpoints=sorted(ids);pair=digest(endpoints)
                if pair in pairs and pairs[pair]!=endpoints:raise ValueError('Pair hash collision')
                pairs[pair]=endpoints
                writer.writerow(dict(row,boundary=boundary,interval_left=ids[0],interval_right=ids[1],domain_pair_key=pair));counts[row['policy']+':'+boundary]+=1
    if source_rows!=sum(r['single_copy_domain_matches_by_policy'].values()):raise ValueError('Missing domain match rows')
    with (out/'intervals.jsonl').open('w') as f:
        for key in sorted(intervals):f.write(json.dumps(intervals[key],separators=(',',':'))+'\n')
    with (out/'domain_pairs.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['domain_pair_key','interval_a','interval_b'])
        w.writerows((key,*pairs[key]) for key in sorted(pairs))
    verify();result=dict(status='complete_duplication_domain_pair_inventory_pending_readback',plan_sha256=ph,source_control_receipt_sha256=sha(source/'receipt.json'),source_domain_match_rows=source_rows,policy_boundary_links=sum(counts.values()),link_counts=dict(counts),unique_intervals=len(intervals),unique_interval_pairs=len(pairs),unique_models=len({(v['model_id'],v['version']) for v in intervals.values()}),interval_residues=sum(v['length'] for v in intervals.values()),maximum_interval_length=max(v['length'] for v in intervals.values()),planned_directed_mask_dispositions=4*len(pairs),artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every verified single-copy candidate Domain match expanded to alignment and envelope boundaries. All policy, source protein pair, scope, accession and hit links retained; only identical coordinate-interval computations deduplicated. No domain coordinates extracted or compared yet; residue confidence, PAE, independent output readback and biological interpretation remain pending.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','scope','link_counts']},indent=2))


if __name__=='__main__':main()
