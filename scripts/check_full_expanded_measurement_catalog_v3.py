#!/usr/bin/env python3
"""Directed catalog contracts, nullable failures and rehashed corrupt exports.

Prior native fits, matching and original journals are synthetic contracts;
the tests cover the new normalization/readback software, not real fitting or
full scientific source acceptance. The complete production scope is unchanged.
"""
import argparse
from collections import Counter
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from background_measurement_union_sources import NUMERIC_INTS,NUMERIC_FLOATS,GEOMETRY_INTS,GEOMETRY_FLOATS
from check_full_matched_coverage_cases import write,table,zipped_table,SCREENS
from summarize_primary_usable_orders_v2 import pair_summary
from run_ortholog_pair_guide_comparison import sha


def invoke(cmd, success=True):
    r=subprocess.run(cmd,text=True,capture_output=True)
    if success and r.returncode:raise RuntimeError(r.stderr+r.stdout)
    if not success:assert r.returncode!=0,'Bad catalog accepted'


def setup_initial_contract(root):
    native,diagnostic,geometry,summary,tquality,bquality,union,queue=[root/n for n in ['native','diagnostic','geometry','summary','tquality','bquality','union','queue']]
    for p in [native,diagnostic,geometry,summary,tquality,bquality,union,queue]:p.mkdir()
    pairends={hashlib.sha256(json.dumps(e,separators=(',',':')).encode()).hexdigest():e for e in [[('A',6),('B',10)],[('C',6),('D',6)]]}
    pairlist=sorted(pairends);native_rows=[];numeric=[];shapes=[];summaries=[];background=[];checkpoints=[];quality_rows={'target':[],'background':[]};targetusable=0
    for role in ['target','background']:
        for i,pair in enumerate(pairlist):
            ends=pairends[pair]
            for mask in ['full','plddt70']:
                dispositions={};numbers={};geos={};usability=[]
                q=dict(pair_key=pair,mask=mask,model_a=ends[0][0],version_a=ends[0][1],model_b=ends[1][0],version_b=ends[1][1],length_a=100,length_b=120)
                for order in [0,1]:
                    status='aligned';n=80;rmsd=0.0 if i==1 else 1.5+order;degenerate=False
                    if role=='target' and mask=='plddt70':
                        if i==0 and order==0:degenerate=True;n=10
                        else:status='parse_error' if i==0 else 'input_unavailable'
                    if role=='background':
                        if i==0 and mask=='plddt70':
                            if order==0:status='timeout'
                            else:n=2;degenerate=True
                        elif i==1 and mask=='full':status='native_error' if order==0 else 'input_unavailable'
                        elif i==1:rmsd=2.0
                    number=shape=None;why=[]
                    if status=='aligned':
                        number={k:0 for k in NUMERIC_INTS};number.update({k:0.25 for k in NUMERIC_FLOATS});number.update(aligned_length=n,joint_plddt70_pairs=n,rmsd_recomputed=rmsd,rmsd_native=rmsd,rmsd_rounding_error=0.0,rmsd_status='within_printed_rounding',coverage_left=n/(100 if order==0 else 120),coverage_right=n/(120 if order==0 else 100),joint_plddt70_fraction=1.0)
                        shape={k:1 for k in GEOMETRY_INTS};shape.update({k:1.0 for k in GEOMETRY_FLOATS});shape.update(aligned_length=n,geometry_status='degenerate_at_numeric_tolerance' if degenerate else 'unique_at_numeric_tolerance')
                        if n<3:why.append('fewer_than_three_pairs')
                        if degenerate:why.append('nonunique_rotation')
                    else:why=[status]
                    usable=not why;usability.append(usable);dispositions[order]=status
                    source_order=1-order if role=='background' and i==0 else order
                    cp=(native/'pairs'/pair[:2] if role=='target' else union)/f'{pair}-{mask}-{source_order}.json';write(cp,dict(synthetic_native_contract=True,role=role,pair_key=pair,mask=mask,order=source_order,status=status));checkpoints.append(cp)
                    if role=='target':
                        native_rows.append(dict(path=str(cp.relative_to(native)),sha256=sha(cp),status=status));targetusable+=usable
                        if number is not None:
                            nr=dict(pair_key=pair,mask=mask,order=order,**{k:str(v) for k,v in number.items()});gr=dict(pair_key=pair,mask=mask,order=order,rmsd_status=number['rmsd_status'],**{k:str(v) for k,v in shape.items()});numeric.append(nr);shapes.append(gr);numbers[order]=nr;geos[order]=gr
                        q.update({f'order{order}_status':'aligned' if usable else 'excluded_numerically' if status=='aligned' else status,f'order{order}_native_status':status,f'order{order}_numerical_exclusion_reasons':';'.join(why) if status=='aligned' else '',f'order{order}_aligned_length':str(float(n)) if usable else ''})
                    else:
                        background.append(dict(pair_key=pair,mask=mask,order=order,directed_endpoints=[dict(model_id=e[0],version=e[1]) for e in (ends if order==0 else ends[::-1])],selected_source='reference_old' if i==0 else 'new_native',source_checkpoint=str(cp),source_checkpoint_sha256=sha(cp),source_order=source_order,source_native_status=status,numerical=number,geometry=shape,numerical_usable=usable,numerical_exclusion_reasons=why))
                        q.update({f'order{order}_native_status':status,f'order{order}_numerical_usable':int(usable),f'order{order}_aligned_length':n if number is not None else ''})
                if role=='target':summaries.append(pair_summary(pair,mask,dispositions,numbers,geos))
                for spec in SCREENS:
                    passed=all(usability) and spec['minimum_aligned_residues']<=80 and spec['minimum_original_coverage']<=80/120
                    q[spec['id']+'_pass']=int(passed);q[spec['id']+'_exclusions']='' if passed else 'synthetic_screen_exclusion'
                quality_rows[role].append(q)
    for p,rows in [(native/'checkpoint_manifest.tsv',native_rows),(diagnostic/'numeric_readback.tsv',numeric),(geometry/'alignment_geometry.tsv',shapes),(summary/'pair_mask_order_summary.tsv',summaries),(tquality/'pair_mask_coverage.tsv',quality_rows['target']),(bquality/'pair_mask_coverage.tsv',quality_rows['background'])]:table(p,rows)
    table(queue/'model_pairs.tsv',[dict(pair_key=p,model_a=e[0][0],version_a=e[0][1],model_b=e[1][0],version_b=e[1][1]) for p,e in pairends.items()])
    states=union/'background_measurement_dispositions.jsonl.gz'
    with gzip.open(states,'wt') as h:
        for row in background:h.write(json.dumps(row)+'\n')
    np,sp,sa=root/'native-plan.json',root/'summary-plan.json',root/'summary-readback.json'
    write(np,dict(output=str(native),queue=str(queue)));write(sp,dict(output=str(summary),native=str(native),diagnostic=str(diagnostic),geometry=str(geometry)))
    write(summary/'receipt.json',dict(status='complete_primary_numerically_usable_order_summary_pending_independent_readback',plan_sha256=sha(sp),pairs=2,pair_mask_rows=4));write(sa,dict(status='passed_full_primary_usable_order_summary_readback',source_receipt_sha256=sha(summary/'receipt.json'),pair_mask_rows=4,numerically_usable_directions=targetusable))
    up,uc=root/'union-plan.json',root/'union-completed.json';write(up,dict(output=str(union)));write(union/'receipt.json',dict(status='complete_full_background_measurement_union_pending_independent_readback',directed_dispositions=8))
    write(uc,dict(status='complete_verified_full_background_measurement_union',full_pairs=2,producer_receipt=str(union/'receipt.json'),producer_receipt_sha256=sha(union/'receipt.json')))
    tp,bp=root/'target-plan.json',root/'background-plan.json';write(tp,dict(output=str(tquality),summary=str(summary),summary_audit=str(sa),screens=SCREENS));write(bp,dict(output=str(bquality),measurement_plan=str(up),measurement_completion=str(uc),screens=SCREENS))
    mp=root/'matching-plan.json';write(mp,dict(target_plan=str(tp),background_plan=str(bp),screens=SCREENS))
    bindings={str(p):sha(p) for p in [np,sp,sa,up,uc,tp,bp,mp,summary/'receipt.json',union/'receipt.json',native/'checkpoint_manifest.tsv',queue/'model_pairs.tsv',diagnostic/'numeric_readback.tsv',geometry/'alignment_geometry.tsv',summary/'pair_mask_order_summary.tsv',tquality/'pair_mask_coverage.tsv',bquality/'pair_mask_coverage.tsv',states,*checkpoints]}
    archive=root/'matching-archive.json';write(archive,dict(status='complete_verified_full_matched_coverage_archive',summary=dict(selected_records=7),services=[dict(synthetic_fixture_only=True)]*2,source_hashes=bindings))
    mc=root/'matching-completed.json';write(mc,dict(status='complete_verified_full_fixed_matched_coverage_attrition',selected_records=7,scientific_eligibility=False,source_plan=str(mp),source_plan_sha256=sha(mp),exact_process_journals_checked=2,full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(bindings)))
    pp=root/'catalog-plan.json';write(pp,dict(matching_completion=str(mc),matching_plan=str(mp),target_summary_plan=str(sp),target_summary_readback=str(sa),target_native_plan=str(np),expected=dict(target_pairs=2,background_pairs=2,selected_records=7,target_numeric_states=len(numeric)),screens=SCREENS,pins={},resources=dict(minimum_free_disk_gib=0),output=str(root/'catalog'),scope=__doc__));return pp


def setup(root):
    pp=setup_initial_contract(root);plan=json.loads(pp.read_text())
    sp=json.loads(Path(plan['target_summary_plan']).read_text());native,diagnostic,geometry,summary=[Path(sp[k]) for k in ['native','diagnostic','geometry','output']]
    manifest=list(csv.DictReader((native/'checkpoint_manifest.tsv').open(),delimiter='\t'))
    counts=dict(Counter(Path(r['path']).stem.rsplit('-',2)[1]+':'+r['status'] for r in manifest))
    write(native/'receipt.json',dict(status='complete_duplication_alignment_dispositions_pending_readback',plan_sha256=sha(plan['target_native_plan']),distinct_model_pairs=2,directed_dispositions=8,counts=counts,artifacts={'checkpoint_manifest.tsv':sha(native/'checkpoint_manifest.tsv')}))
    write(diagnostic/'receipt.json',dict(status='complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance',producer_receipt_sha256=sha(native/'receipt.json'),directed_dispositions=8,numerically_checked_alignments=plan['expected']['target_numeric_states'],counts=counts,artifacts={'numeric_readback.tsv':sha(diagnostic/'numeric_readback.tsv')}))
    write(geometry/'receipt.json',dict(status='complete_primary_diagnostic_geometry_pending_independent_readback',diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),alignments=plan['expected']['target_numeric_states'],artifacts={'alignment_geometry.tsv':sha(geometry/'alignment_geometry.tsv')}))
    sr=json.loads((summary/'receipt.json').read_text());sr.update(native_receipt_sha256=sha(native/'receipt.json'),diagnostic_receipt_sha256=sha(diagnostic/'receipt.json'),geometry_receipt_sha256=sha(geometry/'receipt.json'));write(summary/'receipt.json',sr)
    sa=Path(plan['target_summary_readback']);ar=json.loads(sa.read_text());ar['source_receipt_sha256']=sha(summary/'receipt.json');write(sa,ar)
    mc=Path(plan['matching_completion']);c=json.loads(mc.read_text());archive=Path(c['full_hash_archive']);a=json.loads(archive.read_text())
    for p in [summary/'receipt.json',sa]:a['source_hashes'][str(p)]=sha(p)
    # Deliberately omit the native manifest/receipt from the downstream archive,
    # matching production: the unchanged closed summary binds the original
    # receipt -> manifest -> every native checkpoint chain.
    a['source_hashes'].pop(str(native/'checkpoint_manifest.tsv'))
    write(archive,a);c.update(full_hash_archive_sha256=sha(archive),bound_source_hashes=len(a['source_hashes']));write(mc,c)
    return pp


def run():
    with tempfile.TemporaryDirectory(prefix='expanded-measurement-catalog-fixture-') as temp:
        root=Path(temp);pp=setup(root);out=root/'catalog';producer=[sys.executable,'scripts/export_full_expanded_measurement_catalog_v3.py','--plan',str(pp)];reader=[sys.executable,'scripts/readback_full_expanded_measurement_catalog_v3.py','--plan',str(pp),'--output']
        invoke(producer);invoke(reader+[str(root/'passed.json')]);rp=out/'receipt.json';receipt=json.loads(rp.read_text());assert receipt['directed_states']==16
        rows=list(csv.DictReader(gzip.open(out/'directed_measurements.tsv.gz','rt'),delimiter='\t'))
        assert any(r['native_status']!='aligned' and r['rmsd_recomputed']=='' for r in rows)
        assert any(r['numerical_usable']=='0' and r['rmsd_recomputed']!='' for r in rows)
        assert any(r['numerical_usable']=='1' and r['rmsd_recomputed']=='0.0' for r in rows)
        assert any(r['source_order']!=r['order'] for r in rows)
        invoke(producer,False);original=rp.read_bytes();blob=(out/'directed_measurements.tsv.gz').read_bytes()
        labels=['missing_state','duplicate_state','changed_role','changed_pair','changed_mask','changed_order','changed_model_version','changed_source_order','changed_source_checkpoint','changed_checkpoint_hash','changed_rmsd','changed_native_tm','changed_geometry','cleared_exclusion','promoted_numeric_failure','missing_as_zero','masked_denominator_substitution','changed_coverage_bits','changed_intersection_bits','changed_summary']
        for label in labels:
            modified=copy.deepcopy(rows);r=copy.deepcopy(receipt)
            if label=='missing_state':modified.pop()
            elif label=='duplicate_state':modified.append(copy.deepcopy(modified[0]))
            elif label=='changed_role':modified[0]['role']='background'
            elif label=='changed_pair':modified[0]['pair_key']='0'*64
            elif label=='changed_mask':modified[0]['mask']='plddt70'
            elif label=='changed_order':modified[0]['order']='1'
            elif label=='changed_model_version':modified[0]['version_b']='6'
            elif label=='changed_source_order':modified[0]['source_order']='1'
            elif label=='changed_source_checkpoint':modified[0]['source_checkpoint']='other'
            elif label=='changed_checkpoint_hash':modified[0]['source_checkpoint_sha256']='0'*64
            elif label=='changed_rmsd':modified[0]['rmsd_recomputed']='999'
            elif label=='changed_native_tm':modified[0]['tm_left_native']='999'
            elif label=='changed_geometry':modified[0]['relative_rotation_curvature']='999'
            elif label=='cleared_exclusion':next(x for x in modified if x['numerical_exclusion_reasons'])['numerical_exclusion_reasons']=''
            elif label=='promoted_numeric_failure':next(x for x in modified if x['native_status']=='aligned' and x['numerical_usable']=='0')['numerical_usable']='1'
            elif label=='missing_as_zero':next(x for x in modified if x['native_status']!='aligned')['rmsd_recomputed']='0.0'
            elif label=='masked_denominator_substitution':modified[0]['original_coverage_a']='1.0'
            elif label=='changed_coverage_bits':modified[0]['pair_mask_pass_bits']='63'
            elif label=='changed_intersection_bits':modified[0]['both_masks_pass_bits']='63'
            else:r['usable_states']+=1
            zipped_table(out/'directed_measurements.tsv.gz',modified);r['artifacts']['directed_measurements.tsv.gz']=sha(out/'directed_measurements.tsv.gz');write(rp,r);invoke(reader+[str(root/(label+'.json'))],False)
        rp.write_bytes(original);(out/'directed_measurements.tsv.gz').write_bytes(blob)
        rp.unlink();(out/'directed_measurements.tsv.gz.tmp').write_text('partial');invoke(producer);invoke(reader+[str(root/'recovered.json')])
        assert list(csv.DictReader(gzip.open(out/'directed_measurements.tsv.gz','rt'),delimiter='\t'))==rows
        return dict(status='passed_full_expanded_measurement_catalog_software_contracts',directed_states=16,raw_numeric_failures_quarantined=True,native_errors_and_null_metrics_retained=True,real_numeric_zero_retained=True,source_order_reversal_preserved=True,shared_physical_pairs_across_roles_not_collapsed=True,full_length_denominators_and_masks_checked=True,completed_restart_refused=True,interrupted_full_replay_passed=True,rejected_rehashed_exports=labels,synthetic_prior_native_matching_and_journal_contracts=True,production_fits_or_journals_tested=False,scope=__doc__)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run();r['source_hashes']={n:sha(n) for n in ['scripts/full_expanded_measurement_catalog_sources_v3.py','scripts/export_full_expanded_measurement_catalog_v3.py','scripts/readback_full_expanded_measurement_catalog_v3.py',__file__]}
    with a.output.open('x') as h:h.write(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
