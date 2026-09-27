#!/usr/bin/env python3
"""Reconstruct the complete architecture-support grid with independent vector counts."""
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from run_ortholog_pair_guide_comparison import sha

POLICIES=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']


def classify(a,b):
    first,second=a[0],b[0]
    empty=[len(first)==0,len(second)==0]
    if all(empty):return 'neither_annotated',None
    if any(empty):return 'one_unannotated',None
    if first!=second:return 'different_ordered_annotations',None
    if [a[1],b[1]]!=[True,True]:return 'shared_but_nonconservative',None
    return 'shared_conservative_architecture',first


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);a=ap.parse_args();config=json.loads(a.config.read_text());ch=sha(a.config)
    plan=json.loads(Path(config['source_plan']).read_text());root=Path(plan['output']);receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['status']=='complete_architecture_matched_support_pending_readback' and receipt['plan_sha256']==sha(config['source_plan'])
    bindings={**plan['pins'],**config['pins'],str(a.config):ch,str(root/'receipt.json'):sha(root/'receipt.json')}
    bindings.update({str(root/name):h for name,h in receipt['artifacts'].items()})
    def verify():
        for p,h in bindings.items():assert sha(p)==h,p
    verify()
    annotations={}
    for path in plan['annotation_tables']:
        with Path(path).open() as f:
            for line in f:
                row=json.loads(line);key=(row['model_id'],row['version'])
                item={policy:(tuple((h['pfam_accession'],h['pfam_type']) for h in hits),len(hits)>0 and {h['conservative_architecture'] for h in hits}=={1}) for policy,hits in row['policies'].items()}
                assert key not in annotations or annotations[key]==item;annotations[key]=item
    targets={}
    with Path(plan['target_links']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=(row['guide'],row['gene_a'],row['gene_b']);assert key not in targets;targets[key]=row
    groups=defaultdict(list);background_counts=Counter()
    with Path(plan['background_links']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            flags=[row[g+'_native_ortholog']=='1' and row[g+'_candidate_status']=='cross_taxon_unreported_candidate' for g in ['profile','mafft']]
            for i,guide in enumerate(['profile','mafft']):
                included={'guide_native_ortholog':flags[i],'both_guides_native_ortholog':all(flags),'both_guides_unreported_parents':all(flags) and row['both_parents_unreported']=='1'}
                for setting,keep in included.items():
                    if not keep:continue
                    left=annotations[row['model_id_a'],int(row['version_a'])];right=annotations[row['model_id_b'],int(row['version_b'])]
                    for policy in POLICIES:
                        status,sig=classify(left[policy],right[policy]);background_counts['|'.join([guide,setting,policy,status])]+=1
                        if sig is not None:groups[guide,setting,policy,row[guide+'_family'],sig].append((float(row[guide+'_sequence_pair_distance']),row['taxon_a'],row['taxon_b']))
    arrays={k:(np.array([r[0] for r in rows]),np.array([r[1] for r in rows]),np.array([r[2] for r in rows])) for k,rows in groups.items()}
    empty=(np.array([],dtype=float),np.array([],dtype=str),np.array([],dtype=str));counts=Counter();source_seen=set();n=0
    with Path(plan['support_table']).open() as src,(root/'architecture_support.tsv').open() as exported:
        actual=iter(csv.DictReader(exported,delimiter='\t'))
        for row in csv.DictReader(src,delimiter='\t'):
            key=(row['guide'],row['gene_a'],row['gene_b']);target=targets[key];setting=row['background_set'];identity=key+(setting,);assert identity not in source_seen;source_seen.add(identity)
            assert row['family']==target['family'] and row['gene_node']==target['gene_node'] and row['taxon_id']==target['taxon_id']
            left=annotations[target['model_a'],int(target['version_a'])];right=annotations[target['model_b'],int(target['version_b'])]
            for policy in POLICIES:
                status,sig=classify(left[policy],right[policy]);values,ta,tb=arrays.get((row['guide'],setting,policy,row['family'],sig),empty);focal=(ta==row['taxon_id'])|(tb==row['taxon_id']);d=float(row['target_sequence_distance'])
                expected={k:row[k] for k in ['guide','family','taxon_id','gene_a','gene_b','gene_node','target_same_model','target_sequence_distance','background_set','distance_rule']}
                expected.update(policy=policy,target_architecture_status=status,ordered_signature=json.dumps(sig,separators=(',',':')) if sig is not None else '',architecture_backgrounds=str(len(values)),focal_architecture_backgrounds=str(int(focal.sum())),nonfocal_architecture_backgrounds=str(int((~focal).sum())))
                prefix='|'.join([row['guide'],setting,policy]);counts[prefix+'|'+status]+=1;counts[prefix+'|targets']+=1;counts[prefix+'|targets_with_architecture_background']+=int(len(values)>0)
                for factor,label in [(1.25,'1_25'),(1.5,'1_5'),(2.,'2_0')]:
                    mask=(values>=d/factor)&(values<=d*factor);total=int(mask.sum());local=int((mask&focal).sum())
                    expected['within_factor_'+label]=str(total);expected['focal_within_factor_'+label]=str(local);expected['nonfocal_within_factor_'+label]=str(total-local)
                    counts[prefix+'|targets_within_factor_'+label]+=int(total>0);counts[prefix+'|targets_focal_within_factor_'+label]+=int(local>0)
                assert next(actual,None)==expected,(key,setting,policy);n+=1
            if n%200000==0:print(n,'rows verified',flush=True)
        assert next(actual,None) is None
    assert len(source_seen)==len(targets)*3 and n==len(targets)*12==receipt['rows'] and len(targets)==receipt['targets']
    assert dict(counts)==receipt['counts'] and dict(background_counts)==receipt['background_counts'];verify()
    result=dict(status='passed_full_architecture_matched_support_readback',config_sha256=ch,producer_receipt_sha256=sha(root/'receipt.json'),targets=len(targets),rows=n,checker_sha256=sha(__file__),scope='Every target/set/policy row, exact signature and status, full background qualification pool and family/focal/nonfocal distance count independently reconstructed. Unsorted NumPy masks replace producer binary interval counts; missing/different/nonconservative architecture targets retained. This verifies support, not final matching, structural responses or biological effects.')
    p=Path(config['output']);assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
