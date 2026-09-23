#!/usr/bin/env python3
"""Compare paired domain boundaries within one completed structural partition."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import shutil
import time
import psutil
from catalog_whole_proteome_structures import sha


def disposition(alignment, envelope, assignments):
    a=assignments.get(alignment);e=assignments.get(envelope)
    if a is None and e is None:return 'both_unclustered'
    if a is None:return 'alignment_unclustered'
    if e is None:return 'envelope_unclustered'
    if alignment==envelope:return 'identical_interval'
    return 'distinct_intervals_same_cluster' if a==e else 'distinct_intervals_different_clusters'


def compare(links, assignments):
    pairs={}
    for row in links:
        key=row['model_key'],row['hit_id'];bounds=pairs.setdefault(key,{})
        if row['boundary'] not in {'alignment','envelope'} or row['boundary'] in bounds:raise ValueError('Wrong/repeated boundary')
        bounds[row['boundary']]=row['interval_id']
    rows=[]
    for (model,hit),bounds in sorted(pairs.items()):
        if set(bounds)!={'alignment','envelope'}:raise ValueError('Incomplete paired boundaries')
        a,e=bounds['alignment'],bounds['envelope']
        rows.append(dict(model_key=model,hit_id=hit,alignment_interval=a,envelope_interval=e,
                         alignment_cluster=assignments.get(a,''),envelope_cluster=assignments.get(e,''),
                         disposition=disposition(a,e,assignments)))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed input: '+p)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    def state(s): (out/'state.json').write_text(json.dumps({'status':s})+'\n')
    state('waiting_for_completed_domain_clustering')
    while True:
        try:
            p=psutil.Process(plan['predecessor']['pid']);live=p.create_time()==plan['predecessor']['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify()
    if shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient free disk')
    clusters=Path(plan['clusters']);crp=clusters/'receipt.json';cr=json.loads(crp.read_text())
    db=Path(plan['database']);drp=db/'receipt.json';dr=json.loads(drp.read_text())
    if cr['status']!='complete_domain_candidate_partition_membership_readback' or cr['plan_sha256']!=sha(plan['cluster_plan']) or cr['database_receipt_sha256']!=sha(drp):raise ValueError('Incomplete/unbound clustering')
    if dr['status']!='complete_domain_foldseek_database_with_full_sequence_and_coordinate_readback' or dr['plan_sha256']!=sha(plan['database_plan']):raise ValueError('Incomplete/unbound database')
    lookup=db/'domains.lookup';members=clusters/'cluster_members.tsv'
    for path,digest in [(lookup,dr['artifacts'][lookup.name]),(members,cr['artifacts'][members.name])]:
        if sha(path)!=digest:raise ValueError('Changed partition input')
    expected=set()
    with lookup.open() as f:
        for line in f:
            _,name,_=line.rstrip('\n').split('\t')
            if name in expected:raise ValueError('Repeated lookup identity')
            expected.add(name)
    assignments={};representatives=set()
    with members.open() as f:
        for line in f:
            rep,member=line.rstrip('\n').split('\t')
            if member in assignments or rep not in expected or member not in expected:raise ValueError('Invalid partition identity')
            assignments[member]=rep;representatives.add(rep)
    if set(assignments)!=expected or any(assignments[r]!=r for r in representatives) or len(expected)!=cr['intervals'] or len(expected)!=dr['models']:raise ValueError('Incomplete partition')
    manifest=Path(plan['manifest']);mr=json.loads((manifest/'receipt.json').read_text())
    if mr['status']!='complete_all_candidate_domain_extraction_manifest':raise ValueError('Incomplete manifest')
    with gzip.open(manifest/'boundary_links.tsv.gz','rt') as f:links=list(csv.DictReader(f,delimiter='\t'))
    universe={r['interval_id'] for r in links}
    if len(links)!=mr['boundary_links'] or len(universe)!=mr['unique_intervals'] or not expected<=universe or len(universe-expected)!=cr['excluded_rejected_intervals']:raise ValueError('Export/manifest scope differs')
    rows=compare(links,assignments)
    if len(rows)!=mr['candidate_model_hit_pairs']:raise ValueError('Wrong candidate pair count')
    target=out/'boundary_cluster_dispositions.tsv.gz'
    with gzip.open(target,'wt') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    counts=dict(Counter(r['disposition'] for r in rows));verify()
    if sha(lookup)!=dr['artifacts'][lookup.name] or sha(members)!=cr['artifacts'][members.name]:raise ValueError('Partition changed during comparison')
    result={'status':'complete_all_domain_boundary_cluster_dispositions','candidate_pairs':len(rows),'intervals':len(universe),'clustered_intervals':len(expected),'unclustered_intervals':len(universe-expected),'dispositions':counts,'plan_sha256':ph,'cluster_receipt_sha256':sha(crp),'database_receipt_sha256':sha(drp),'artifacts':{target.name:sha(target)},'scope':'Every original model/hit boundary pair compared within one candidate partition, retaining excluded intervals. Identical intervals are separated from distinct boundaries sharing a cluster. Disagreement is sensitivity within this partition, not separate-run stability, direct structural distance, independent transitions, homology or biological divergence. Full output readback and parameter sensitivity remain required.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');state('complete')


if __name__=='__main__':main()
