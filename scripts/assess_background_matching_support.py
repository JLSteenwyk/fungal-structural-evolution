#!/usr/bin/env python3
"""Measure pre-outcome family/taxon/sequence-distance support for every duplicate target."""
import argparse,bisect,csv,json,math,time
from collections import Counter,defaultdict
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha

FACTORS=(1.25,1.5,2.0)
SETS=('guide_native_ortholog','both_guides_native_ortholog','both_guides_unreported_parents')


def interval_count(values,distance,factor):
    assert math.isfinite(distance) and distance>=0 and factor>=1
    low=distance/factor;high=distance*factor
    return bisect.bisect_right(values,high)-bisect.bisect_left(values,low)


def supported_sets(row,guide):
    local=row[guide+'_candidate_status']=='cross_taxon_unreported_candidate' and row[guide+'_native_ortholog']=='1'
    if not local:return []
    result=[SETS[0]]
    if row['modeled_candidate_in_both_guides']=='1' and row['profile_native_ortholog']==row['mafft_native_ortholog']=='1':
        result.append(SETS[1])
        if row['both_parents_unreported']=='1':result.append(SETS[2])
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for p,h in plan['pins'].items():assert sha(p)==h,p
    verify();dep=plan['dependency']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();root=Path(plan['source']);r=json.loads((root/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text());rh=sha(root/'receipt.json');ah=sha(plan['readback'])
    assert audit['status']=='passed_full_background_orthology_membership_readback' and audit['producer_receipt_sha256']==rh
    table=root/'candidate_orthology_membership.tsv';assert sha(table)==r['artifacts'][table.name]
    with table.open() as f:background=list(csv.DictReader(f,delimiter='\t'))
    assert len(background)==audit['candidate_rows']
    inv=Path(plan['inventory']);ir=json.loads((inv/'receipt.json').read_text());ia=json.loads(Path(plan['inventory_readback']).read_text())
    assert ia['status']=='passed_full_terminal_sister_inventory_readback' and ia['producer_receipt_sha256']==sha(inv/'receipt.json')
    review=Path(plan['targets']);tr=json.loads((review/'receipt.json').read_text());assert tr['status']=='complete_duplication_candidate_join_and_tree_review'
    out=Path(plan['output']);out.mkdir(exist_ok=False);summaries=[];total=0
    with (out/'target_support.tsv').open('w') as handle:
        writer=None
        for guide in ['profile','mafft']:
            pools=defaultdict(list);focal=defaultdict(list);model_classes=Counter()
            for row in background:
                for setting in supported_sets(row,guide):
                    family=row[guide+'_family'];distance=float(row[guide+'_sequence_pair_distance']);assert math.isfinite(distance) and distance>=0
                    assert row['taxon_a']!=row['taxon_b']
                    pools[setting,family].append(distance)
                    for taxon in [row['taxon_a'],row['taxon_b']]:focal[setting,family,taxon].append(distance)
                    model_classes[setting+':'+row[guide+'_model_coverage']]+=1
            for values in list(pools.values())+list(focal.values()):values.sort()
            path=review/(guide+'_candidate_tree_checks.tsv');assert sha(path)==tr['artifacts'][path.name]
            with path.open() as f:targets={(row['gene_a'],row['gene_b']):row for row in csv.DictReader(f,delimiter='\t')}
            assert len(targets)==sum(1 for _ in path.open())-1
            seen=set();counts=Counter();p=inv/(guide+'_terminal_sisters.tsv');assert sha(p)==ir['artifacts'][p.name]
            with p.open() as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    key=(row['gene_a'],row['gene_b'])
                    if key not in targets:continue
                    assert key not in seen;seen.add(key);target=targets[key]
                    assert target['tree_status']=='exact_reported_pair' and row['candidate_status']=='reported_duplication'
                    assert row['taxon_a']==row['taxon_b']==target['taxon_id'] and row['family']==target['family'] and row['gene_node']==target['gene_node']
                    d=float(row['sequence_pair_distance']);taxon=target['taxon_id'];family=target['family']
                    for setting in SETS:
                        values=pools[setting,family];local=focal[setting,family,taxon]
                        record=dict(guide=guide,family=family,taxon_id=taxon,gene_a=key[0],gene_b=key[1],gene_node=row['gene_node'],target_same_model=target['same_model'],target_sequence_distance=d,background_set=setting,distance_rule='zero_exact_only' if d==0 else 'multiplicative_distance_range',family_backgrounds=len(values),focal_taxon_backgrounds=len(local),nonfocal_backgrounds=len(values)-len(local))
                        for factor in FACTORS:
                            label=str(factor).replace('.','_');n=interval_count(values,d,factor);nf=interval_count(local,d,factor)
                            record['within_factor_'+label]=n;record['focal_within_factor_'+label]=nf;record['nonfocal_within_factor_'+label]=n-nf
                            counts[setting+':targets_within_factor_'+label]+=int(n>0)
                            counts[setting+':targets_focal_within_factor_'+label]+=int(nf>0)
                        counts[setting+':targets']+=1;counts[setting+':no_family_background']+=int(not values);counts[setting+':zero_target_distance']+=int(d==0)
                        if writer is None:writer=csv.DictWriter(handle,list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
                        writer.writerow(record);total+=1
            assert seen==set(targets)
            summaries.append(dict(guide=guide,targets=len(targets),counts=dict(counts),background_model_classes=dict(model_classes)))
            print(json.dumps(summaries[-1]),flush=True)
    verify();assert sha(root/'receipt.json')==rh and sha(plan['readback'])==ah and sha(table)==r['artifacts'][table.name]
    result=dict(status='complete_background_matching_support_pending_readback',plan_sha256=ph,source_receipt_sha256=rh,source_readback_sha256=ah,rows=total,guides=summaries,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every reviewed two-model duplicate target, three background qualification sets, same-family pool counts and focal/nonfocal taxon overlap. Sequence-distance factors 1.25/1.5/2 are descriptive sensitivity ranges, zero targets require exact-zero backgrounds. Identical models and no-support targets retained. No matches selected; domain/confidence/length/phylogenetic balance and structural responses not yet assessed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
