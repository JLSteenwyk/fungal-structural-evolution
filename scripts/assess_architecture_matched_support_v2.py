#!/usr/bin/env python3
"""Assess same-family background support within conservatively shared architectures."""
import argparse,bisect,csv,json
from collections import Counter,defaultdict
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

POLICIES=('alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore')


def pair_signature(left,right):
    a,b=left[0],right[0]
    if not a and not b:return 'neither_annotated',None
    if not a or not b:return 'one_unannotated',None
    if a!=b:return 'different_ordered_annotations',None
    if not left[1] or not right[1]:return 'shared_but_nonconservative',None
    return 'shared_conservative_architecture',a


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    source_bindings={}
    def verify():
        assert sha(a.plan)==ph
        for p,h in {**plan['pins'],**source_bindings}.items():assert sha(p)==h,p
    verify()
    for item in plan['verified_sources']:
        r=json.loads(Path(item['receipt']).read_text());audit=json.loads(Path(item['readback']).read_text())
        assert audit['status']==item['readback_status'] and audit['producer_receipt_sha256']==sha(item['receipt'])
        source_bindings[item['receipt']]=sha(item['receipt'])
        source_bindings[item['readback']]=sha(item['readback'])
        for name,h in r['artifacts'].items():
            path=str(Path(item['receipt']).parent/name)
            assert sha(path)==h,path
            if path in source_bindings:assert source_bindings[path]==h
            source_bindings[path]=h
    required=[*plan['annotation_tables'],plan['background_links'],plan['target_links'],plan['support_table']]
    for path in required:assert path in source_bindings or path in plan['pins'],('Unbound matching input',path)
    verify()
    annotations={}
    for path in plan['annotation_tables']:
        with Path(path).open() as f:
            for line in f:
                row=json.loads(line);key=(row['model_id'],row['version']);values={}
                for policy,hits in row['policies'].items():
                    signature=tuple((h['pfam_accession'],h['pfam_type']) for h in hits)
                    values[policy]=(signature,bool(hits) and all(h['conservative_architecture']==1 for h in hits))
                if key in annotations:assert annotations[key]==values
                annotations[key]=values
    targets={}
    with Path(plan['target_links']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=(row['guide'],row['gene_a'],row['gene_b']);assert key not in targets
            targets[key]=row
    pools=defaultdict(list);focal=defaultdict(list);background_counts=Counter()
    with Path(plan['background_links']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            left=annotations[row['model_id_a'],int(row['version_a'])];right=annotations[row['model_id_b'],int(row['version_b'])]
            for guide in ['profile','mafft']:
                if row[guide+'_candidate_and_native_ortholog']!='1':continue
                settings=['guide_native_ortholog']
                if row['candidate_and_native_ortholog_both']=='1':settings.append('both_guides_native_ortholog')
                if row['both_guides_and_parents_unreported']=='1':settings.append('both_guides_unreported_parents')
                family=row[guide+'_family'];d=float(row[guide+'_sequence_pair_distance'])
                for policy in POLICIES:
                    status,sig=pair_signature(left[policy],right[policy])
                    for setting in settings:
                        background_counts[guide+'|'+setting+'|'+policy+'|'+status]+=1
                        if sig is None:continue
                        key=(guide,setting,policy,family,sig);pools[key].append(d)
                        for taxon in [row['taxon_a'],row['taxon_b']]:focal[key+(taxon,)].append(d)
    for values in list(pools.values())+list(focal.values()):values.sort()
    out=Path(plan['output']);out.mkdir(exist_ok=False);counts=Counter();seen=set();total=0
    with Path(plan['support_table']).open() as src,(out/'architecture_support.tsv').open('w') as dest:
        writer=None
        for row in csv.DictReader(src,delimiter='\t'):
            key=(row['guide'],row['gene_a'],row['gene_b']);target=targets[key]
            assert row['family']==target['family'] and row['taxon_id']==target['taxon_id'] and row['gene_node']==target['gene_node']
            record_key=key+(row['background_set'],);assert record_key not in seen;seen.add(record_key)
            left=annotations[target['model_a'],int(target['version_a'])];right=annotations[target['model_b'],int(target['version_b'])]
            for policy in POLICIES:
                status,sig=pair_signature(left[policy],right[policy]);index=(row['guide'],row['background_set'],policy,row['family'],sig)
                values=pools.get(index,[]);local=focal.get(index+(row['taxon_id'],),[]);d=float(row['target_sequence_distance'])
                record={k:row[k] for k in ['guide','family','taxon_id','gene_a','gene_b','gene_node','target_same_model','target_sequence_distance','background_set','distance_rule']}
                record.update(policy=policy,target_architecture_status=status,ordered_signature=json.dumps(sig,separators=(',',':')) if sig is not None else '',architecture_backgrounds=len(values),focal_architecture_backgrounds=len(local),nonfocal_architecture_backgrounds=len(values)-len(local))
                prefix=row['guide']+'|'+row['background_set']+'|'+policy
                counts[prefix+'|'+status]+=1;counts[prefix+'|targets']+=1;counts[prefix+'|targets_with_architecture_background']+=int(bool(values))
                for factor,label in [(1.25,'1_25'),(1.5,'1_5'),(2.,'2_0')]:
                    n=bisect.bisect_right(values,d*factor)-bisect.bisect_left(values,d/factor);nf=bisect.bisect_right(local,d*factor)-bisect.bisect_left(local,d/factor)
                    assert n<=int(row['within_factor_'+label]) and nf<=int(row['focal_within_factor_'+label])
                    record['within_factor_'+label]=n;record['focal_within_factor_'+label]=nf;record['nonfocal_within_factor_'+label]=n-nf
                    counts[prefix+'|targets_within_factor_'+label]+=int(n>0);counts[prefix+'|targets_focal_within_factor_'+label]+=int(nf>0)
                if writer is None:writer=csv.DictWriter(dest,list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
                writer.writerow(record);total+=1
            if total%200000==0:print(total,'architecture support rows',flush=True)
    assert len(seen)==len(targets)*3 and total==len(targets)*12
    verify();result=dict(status='complete_architecture_matched_support_pending_readback',plan_sha256=ph,targets=len(targets),rows=total,counts=dict(counts),background_counts=dict(background_counts),artifacts={p.name:sha(p) for p in out.iterdir()},source_hashes={**plan['pins'],**source_bindings},scope='Full duplicate target grid, three native-orthology qualification sets and four annotation policies. Background and both target proteins must share the same nonempty conservative ordered annotation signature. Unannotated/different/nonconservative targets retained explicitly, with zero support for this within-architecture estimand. Family, sequence-distance and focal/nonfocal counts only; not final matched controls or causal duplication effects.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','targets','rows']}),flush=True)


if __name__=='__main__':main()
