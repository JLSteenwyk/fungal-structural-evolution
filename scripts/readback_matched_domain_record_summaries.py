#!/usr/bin/env python3
"""Independently rebuild domain-to-record means and all weighting summaries with pandas."""
import argparse,json,hashlib,time,subprocess,itertools
from pathlib import Path
import psutil,numpy as np,pandas as pd

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--launch',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan);dep=json.loads(a.launch.read_text());lh=sha(a.launch);assert dep['plan_sha256']==ph
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',dep['unit'],'-p','ActiveState','-p','ExecMainStatus'],text=True).splitlines());assert state=={'ActiveState':'inactive','ExecMainStatus':'0'}
    assert sha(a.plan)==ph and sha(a.launch)==lh
    for p,h in plan['pins'].items():assert sha(p)==h,p
    root=Path(plan['output']);r=json.loads((root/'receipt.json').read_text());assert r['status']=='complete_matched_domain_record_summaries_pending_independent_readback' and r['plan_sha256']==ph
    for name,h in r['artifacts'].items():assert sha(root/name)==h
    for name in ['measurements','inventory','coverage']:
        assert sha(Path(plan[name])/'receipt.json')==r['source_bindings'][name]['receipt_sha256'] and sha(plan[name+'_audit'])==r['source_bindings'][name]['audit_sha256']
    idcols=['domain_config_id','boundary','mask','target_order','background_order'];flags=[s+'_'+suffix for s in plan['screens'] for suffix in ['mask_pass','both_masks_pass']]
    numeric=['rmsd_target_minus_background','sequence_identity_target_minus_background','target_rmsd_recomputed','background_rmsd_recomputed','target_sequence_identity_exact','background_sequence_identity_exact','target_original_coverage_a','target_original_coverage_b','background_original_coverage_a','background_original_coverage_b','target_aligned_length','background_aligned_length','target_joint_plddt70_fraction','background_joint_plddt70_fraction']
    types={x:str for x in idcols};types.update({x:np.uint8 for x in flags});types.update({x:float for x in numeric})
    d=pd.read_csv(Path(plan['measurements'])/'matched_domain_contrasts.tsv.gz',sep='\t',usecols=idcols+flags+numeric,dtype=types)
    assert len(d)==1200640
    definitions={'rmsd_difference':'rmsd_target_minus_background','identity_difference':'sequence_identity_target_minus_background','target_rmsd':'target_rmsd_recomputed','background_rmsd':'background_rmsd_recomputed','target_identity':'target_sequence_identity_exact','background_identity':'background_sequence_identity_exact'}
    for new,old in definitions.items():d[new]=d[old]
    d['original_coverage_difference']=np.minimum(d.target_original_coverage_a,d.target_original_coverage_b)-np.minimum(d.background_original_coverage_a,d.background_original_coverage_b)
    d['log_aligned_length_ratio']=np.log(d.target_aligned_length)-np.log(d.background_aligned_length)
    d['confidence_fraction_difference']=d.target_joint_plddt70_fraction-d.background_joint_plddt70_fraction
    metrics=list(definitions)+['original_coverage_difference','log_aligned_length_ratio','confidence_fraction_difference'];assert set(metrics)==set(r['metrics'])
    d=d[idcols+flags+metrics]
    selected=pd.read_csv(Path(plan['inventory'])/'selection_domain_links.tsv.gz',sep='\t',usecols=['target_id','background_id','policy','scenario_id','domain_config_id'],dtype=str)
    nodes=[]
    for line in Path(plan['target_nodes']).open():
        x=json.loads(line);nodes.append({k:x[k] for k in ['node_id','guide','family','taxon_id']})
    selected=selected.merge(pd.DataFrame(nodes),left_on='target_id',right_on='node_id',validate='many_to_one');assert len(selected)==2786912
    assert not selected.duplicated(['target_id','policy','scenario_id']).any()
    summary=pd.read_csv(root/'record_summary.tsv',sep='\t',dtype={'target_order':str,'background_order':str});assert len(summary)==82944
    parts=json.loads((root/'partition_manifest.json').read_text());conditions=set();checked=0;config_values=0;groupkeys=['guide','policy','scenario_id']
    scenarios=[x['scenario_id'] for x in json.loads(Path(plan['scenarios']).read_text())];strata=pd.MultiIndex.from_product([['profile','mafft'],plan['policies'],scenarios],names=groupkeys)
    coverage=pd.read_csv(Path(plan['coverage'])/'selected_domain_coverage.tsv',sep='\t').set_index(groupkeys+['boundary','mask_cohort','screen'])
    for part in parts:
        names=['boundary','mask','cohort','screen','target_order','background_order'];condition=tuple(part[k] for k in names);assert condition not in conditions;conditions.add(condition)
        boundary,mask,cohort,screen,to,bo=condition;flag=screen+('_mask_pass' if cohort=='same_mask' else '_both_masks_pass')
        valid=d.boundary.eq(boundary)&d['mask'].eq(mask)&d.target_order.eq(to)&d.background_order.eq(bo)&d[flag].eq(1)
        subset=d.loc[valid,['domain_config_id']+metrics];assert np.isfinite(subset[metrics].to_numpy()).all()
        grouped=subset.groupby('domain_config_id',sort=True);means=grouped[metrics].mean();means['eligible_domains']=grouped.size()
        path=root/part['path'];assert sha(path)==part['sha256'];actual=pd.read_parquet(path).set_index('domain_config_id').sort_index();assert actual.index.is_unique
        assert len(means)==part['configurations'];pd.testing.assert_frame_equal(actual[metrics+['eligible_domains']],means[metrics+['eligible_domains']],check_dtype=False,rtol=1e-10,atol=1e-12);config_values+=len(means)*(len(metrics)+1)
        records=selected.merge(means.reset_index(),on='domain_config_id',validate='many_to_one');g=records.groupby(groupkeys)
        expected=g[metrics].mean().rename(columns={m:m+'_record_mean' for m in metrics})
        for column,label in [('family','family'),('taxon_id','taxon')]:
            result=records.groupby(groupkeys+[column])[metrics].mean().groupby(level=groupkeys).mean().rename(columns={m:m+'_'+label+'_equal_mean' for m in metrics});expected=expected.join(result)
        expected['matched_records']=g.size();expected['families']=g.family.nunique();expected['taxa']=g.taxon_id.nunique();expected['backgrounds']=g.background_id.nunique();expected['mean_eligible_domains']=g.eligible_domains.mean()
        expected=expected.reindex(strata)
        for col in ['matched_records','families','taxa','backgrounds']:expected[col]=expected[col].fillna(0).astype(int)
        chosen=summary.copy()
        for key,value in zip(names,condition):chosen=chosen[chosen[key].eq(value)]
        assert len(chosen)==432;chosen=chosen.set_index(groupkeys).reindex(strata);assert chosen.index.is_unique
        selected_counts=[int(coverage.loc[base+(boundary,mask if cohort=='same_mask' else 'both',screen),'selected_records']) for base in strata];expected['selected_records']=selected_counts
        pd.testing.assert_frame_equal(chosen[expected.columns],expected,check_dtype=False,rtol=1e-10,atol=1e-12);checked+=len(chosen)
        print('Independently summarized',len(conditions),'/ 192',flush=True)
    assert conditions==set(itertools.product(['alignment','envelope'],['full','plddt70'],['same_mask','both_masks'],plan['screens'],['0','1'],['0','1'])) and checked==r['summary_rows']==82944
    for p,h in plan['pins'].items():assert sha(p)==h,p
    result=dict(status='passed_full_matched_domain_record_summary_readback',summary_rows=checked,settings=len(conditions),configuration_values_checked=config_values,source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Every configuration mean rebuilt from full domain measurements and every selected-record/family/taxon mean independently recomputed with pandas, without SQL or producer helpers. Empty strata, counts, denominators and all orders/settings checked. No significance, phylogenetic adjustment or causal inference.')
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
