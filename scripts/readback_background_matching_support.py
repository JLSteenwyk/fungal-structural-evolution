#!/usr/bin/env python3
"""Recount every matching-support row directly with unsorted vector masks."""
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from run_ortholog_pair_guide_comparison import sha

SETS=['guide_native_ortholog','both_guides_native_ortholog','both_guides_unreported_parents']


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_background_matching_support_pending_readback' and receipt['plan_sha256']==sha(a.plan)
    bindings=dict(plan['pins']);bindings[str(a.plan)]=sha(a.plan);bindings[str(root/'receipt.json')]=sha(root/'receipt.json')
    bindings.update({str(root/name):h for name,h in receipt['artifacts'].items()})
    source=Path(plan['source']);sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text())
    assert audit['status']=='passed_full_background_orthology_membership_readback' and audit['producer_receipt_sha256']==receipt['source_receipt_sha256']==sha(source/'receipt.json')
    assert receipt['source_readback_sha256']==sha(plan['readback'])
    for p in [source/'receipt.json',Path(plan['readback'])]:bindings[str(p)]=sha(p)
    table=source/'candidate_orthology_membership.tsv';bindings[str(table)]=sr['artifacts'][table.name]
    def verify():
        for p,h in bindings.items():assert sha(p)==h,p
    verify()
    with table.open() as f:background=list(csv.DictReader(f,delimiter='\t'))
    targets={};pools=defaultdict(list);classes={g:Counter() for g in ['profile','mafft']}
    for guide in ['profile','mafft']:
        for row in background:
            flags=[row[g+'_candidate_status']=='cross_taxon_unreported_candidate' and row[g+'_native_ortholog']=='1' for g in ['profile','mafft']]
            local=flags[['profile','mafft'].index(guide)]
            allowed=[local,all(flags),all(flags) and row['both_parents_unreported']=='1']
            for setting,include in zip(SETS,allowed):
                if include:
                    pools[guide,setting,row[guide+'_family']].append((float(row[guide+'_sequence_pair_distance']),row['taxon_a'],row['taxon_b']))
                    classes[guide][setting+':'+row[guide+'_model_coverage']]+=1
        path=Path(plan['targets'])/(guide+'_candidate_tree_checks.tsv')
        with path.open() as f:review={(r['gene_a'],r['gene_b']):r for r in csv.DictReader(f,delimiter='\t')}
        with (Path(plan['inventory'])/(guide+'_terminal_sisters.tsv')).open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                pair=(row['gene_a'],row['gene_b'])
                if pair not in review:continue
                reviewed=review.pop(pair)
                assert row['family']==reviewed['family'] and row['gene_node']==reviewed['gene_node'] and row['taxon_a']==row['taxon_b']==reviewed['taxon_id']
                assert row['candidate_status']=='reported_duplication'
                targets[(guide,*pair)]=(reviewed,row)
        assert not review
    arrays={key:(np.array([r[0] for r in rows]),np.array([r[1] for r in rows]),np.array([r[2] for r in rows])) for key,rows in pools.items()}
    empty=(np.array([],dtype=float),np.array([],dtype=str),np.array([],dtype=str))
    seen=set();counts={g:Counter() for g in ['profile','mafft']}
    with (root/'target_support.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            guide=row['guide'];key=(guide,row['gene_a'],row['gene_b']);assert key in targets
            setting=row['background_set'];assert setting in SETS and (key,setting) not in seen;seen.add((key,setting))
            reviewed,native=targets[key];d=float(native['sequence_pair_distance']);taxon=reviewed['taxon_id'];family=reviewed['family']
            values,ta,tb=arrays.get((guide,setting,family),empty);focal=(ta==taxon)|(tb==taxon)
            expected=dict(guide=guide,family=family,taxon_id=taxon,gene_a=key[1],gene_b=key[2],gene_node=native['gene_node'],target_same_model=reviewed['same_model'],target_sequence_distance=str(d),background_set=setting,distance_rule='zero_exact_only' if d==0 else 'multiplicative_distance_range',family_backgrounds=str(len(values)),focal_taxon_backgrounds=str(int(focal.sum())),nonfocal_backgrounds=str(int((~focal).sum())))
            for factor,label in [(1.25,'1_25'),(1.5,'1_5'),(2.,'2_0')]:
                mask=(values>=d/factor)&(values<=d*factor);n=int(mask.sum());nf=int((mask&focal).sum())
                expected['within_factor_'+label]=str(n);expected['focal_within_factor_'+label]=str(nf);expected['nonfocal_within_factor_'+label]=str(n-nf)
                counts[guide][setting+':targets_within_factor_'+label]+=int(n>0);counts[guide][setting+':targets_focal_within_factor_'+label]+=int(nf>0)
            assert row==expected,(key,setting)
            counts[guide][setting+':targets']+=1;counts[guide][setting+':no_family_background']+=int(len(values)==0);counts[guide][setting+':zero_target_distance']+=int(d==0)
            if len(seen)%100000==0:print(len(seen),'rows checked',flush=True)
    assert len(seen)==len(targets)*len(SETS)==receipt['rows']
    for summary in receipt['guides']:
        guide=summary['guide'];assert summary['targets']==sum(k[0]==guide for k in targets)
        assert summary['counts']==dict(counts[guide]) and summary['background_model_classes']==dict(classes[guide])
    verify();assert not a.output.exists()
    result=dict(status='passed_full_background_matching_support_readback',producer_receipt_sha256=sha(root/'receipt.json'),plan_sha256=sha(a.plan),target_records=len(targets),support_rows=len(seen),guides=receipt['guides'],checker_sha256=sha(__file__),scope='Every target/set row and all family/focal/nonfocal counts reconstructed from original native target and qualified background rows. Unsorted vector masks independently verify all distance windows, including zero targets; full target grid and summaries checked. Availability is not matched comparability or independent event replication.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','target_records','support_rows']}),flush=True)


if __name__=='__main__':main()
