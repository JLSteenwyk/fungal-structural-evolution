#!/usr/bin/env python3
"""Full 54-scenario descriptive balance/reuse fixture and rehashed export failures.

Prior full matching/physics/journal contracts are synthetic. This runs complete
source I/O, all three baselines, eight features, all masks/screens, node versus
physical-pair reuse, and independent reconstruction. Not a biological pilot.
"""
import argparse
import copy
import csv
import gzip
import itertools
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from check_full_matched_coverage_cases import node, write, table, zipped_table, SCREENS, POLICIES
from full_screened_balance_sources import KEY, MASKS
from full_screened_balance_statistics import FEATURES
from run_ortholog_pair_guide_comparison import sha


def invoke(command, success=True):
    result = subprocess.run(command,text=True,capture_output=True)
    if success and result.returncode: raise RuntimeError(result.stderr+result.stdout)
    if not success: assert result.returncode != 0, 'False screened balance accepted'


def setup(root):
    graph,matching,selection=[root/n for n in ['graph','matching','selection']]
    for p in [graph,matching,selection]:p.mkdir()
    ts=[node('T1','profile','A','B',True),node('T2','profile','A','A',True),node('T3','mafft','A','C',True,True),node('T4','mafft','D','E',True),node('T5','profile','L','M',True),node('T6','profile','N','O',True)]
    bs=[node('B1','profile','F','G'),node('B2','profile','F','H'),node('B3','mafft','I','I'),node('B4','mafft','J','K')]
    b5=copy.deepcopy(bs[0]);b5.update(node_id='B5',gene_a='B5-gene-a',gene_b='B5-gene-b',sequence_distance=0.0);bs.append(b5)
    ts[-2]['sequence_distance']=1.0;ts[-1]['sequence_distance']=0.0
    for i,n in enumerate(ts+bs):
        n.update(mean_ca_plddt_a=70.0+i,mean_ca_plddt_b=60.0+i,fraction_ca_plddt_below50_a=.1,fraction_ca_plddt_below50_b=.2)
    for field in ['mean_ca_plddt_a','mean_ca_plddt_b','fraction_ca_plddt_below50_a','fraction_ca_plddt_below50_b']:b5[field]=bs[0][field]
    target={n['node_id']:n for n in ts};background={n['node_id']:n for n in bs}
    for kind,ns in [('target',ts),('background',bs)]: (graph/(kind+'_nodes.jsonl')).write_text(''.join(json.dumps(n)+'\n' for n in ns))
    tf={'T1':[63,9,9],'T2':[0,0,0],'T3':[63,0,0],'T4':[0,0,0],'T5':[63,63,63],'T6':[63,63,63]};bf={'B1':[63,9,9],'B2':[0,0,0],'B3':[0,0,0],'B4':[63,63,63],'B5':[63,9,9]}
    scenarios=[dict(scenario_id=f'S{i:02d}',synthetic_fixture_only=True) for i in range(1,55)];write(selection/'scenarios.json',scenarios)
    choices=[('T1',POLICIES[0],'S01','B1'),('T5',POLICIES[0],'S01','B1'),('T6',POLICIES[0],'S01','B5'),('T1',POLICIES[0],'S54','B1'),('T2',POLICIES[0],'S01','B2'),('T3',POLICIES[0],'S01','B3'),('T4',POLICIES[0],'S01','B4'),('T1',POLICIES[1],'S01','B2')]
    selected=[]
    for tid,p,sid,bid in choices:
        row=dict(target_id=tid,background_id=bid,guide=target[tid]['guide'],policy=p,scenario_id=sid)
        for j,m in enumerate(MASKS):row.update({f'target_{m}_pass_bits':tf[tid][j],f'control_{m}_pass_bits':bf[bid][j],f'joint_{m}_pass_bits':tf[tid][j]&bf[bid][j]})
        selected.append(row)
    statuses=[]
    for n,p in itertools.product(ts,POLICIES):statuses.append(dict(target_id=n['node_id'],policy=p,guide=n['guide'],target_pair_key=n['pair_key'],target_comparison_disposition='identical_model_no_alignment' if n['same_model'] else 'distinct_model_pair',**{f'target_{m}_pass_bits':tf[n['node_id']][j] for j,m in enumerate(MASKS)}))
    zipped_table(matching/'selected_pair_coverage.tsv.gz',selected);zipped_table(matching/'target_policy_coverage_status.tsv.gz',statuses)
    attrs=[]
    for g,p,s,m,spec in itertools.product(['profile','mafft'],POLICIES,scenarios,MASKS,SCREENS):
        bit=SCREENS.index(spec);mi=MASKS.index(m);group=[r for r in selected if (r['guide'],r['policy'],r['scenario_id'])==(g,p,s['scenario_id'])];allnodes=[n for n in ts if n['guide']==g];cats=Counter((int(bool(tf[r['target_id']][mi]&1<<bit)),int(bool(bf[r['background_id']][mi]&1<<bit))) for r in group);targetpass=sum(bool(tf[n['node_id']][mi]&1<<bit) for n in allnodes);tp=cats[1,0]+cats[1,1]
        row=dict(zip(KEY,[g,p,s['scenario_id'],m,spec['id']]));row.update(all_target_records=len(allnodes),matched_records=len(group),unmatched_records=len(allnodes)-len(group),target_pass_all_records=targetpass,target_pass_matched_records=tp,target_pass_unmatched_records=targetpass-tp,control_pass_matched_records=cats[0,1]+cats[1,1],joint_pass_matched_records=cats[1,1],target_only_pass_matched_records=cats[1,0],control_only_pass_matched_records=cats[0,1],neither_pass_matched_records=cats[0,0]);attrs.append(row)
    table(matching/'matched_attrition_counts.tsv',attrs)
    mp=root/'matching-plan.json';write(mp,dict(output=str(matching),graph=str(graph),selection=str(selection),screens=SCREENS,expected=dict(target_nodes=6,background_nodes=5,selected_records=8,scenarios=54)))
    summary=dict(target_nodes=6,background_nodes=5,target_policy_records=24,scenarios=54,scenario_decisions=1296,selected_records=8,unmatched_decisions=1288,selection_screen_cells=144,full_scenario_screen_cells=23328,attrition_rows=7776,screens=SCREENS,masks=MASKS,guides=['profile','mafft'],policies=POLICIES,node_dispositions={})
    rp,ap=matching/'receipt.json',matching/'readback.json';write(rp,dict(status='complete_full_matched_coverage_pending_independent_readback',plan_sha256=sha(mp),**summary,artifacts={name:sha(matching/name) for name in ['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv']}));write(ap,dict(status='passed_full_matched_coverage_sql_readback',plan_sha256=sha(mp),producer_receipt_sha256=sha(rp),**summary))
    paths=[mp,rp,ap,graph/'target_nodes.jsonl',graph/'background_nodes.jsonl',selection/'scenarios.json']+[matching/name for name in ['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv']];bindings={str(p):sha(p) for p in paths};archive=matching/'archive.json';write(archive,dict(status='complete_verified_full_matched_coverage_archive',summary=summary,services=[dict(synthetic_fixture_only=True)]*2,source_hashes=bindings))
    cp=root/'matching-completed.json';write(cp,dict(status='complete_verified_full_fixed_matched_coverage_attrition',**summary,source_plan=str(mp),source_plan_sha256=sha(mp),producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),independent_readback_sha256=sha(ap),full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(bindings),exact_process_journals_checked=2,scientific_eligibility=False))
    pp=root/'plan.json';write(pp,dict(matching_completion=str(cp),matching_plan=str(mp),screens=SCREENS,guides=['profile','mafft'],policies=POLICIES,features=FEATURES,expected=dict(target_nodes=6,background_nodes=5,selected_records=8,scenarios=54),output=str(root/'balance'),pins={},resources=dict(minimum_free_disk_gib=0),scope=__doc__));return pp


def run():
    with tempfile.TemporaryDirectory(prefix='full-screened-balance-fixture-') as temp:
        root=Path(temp);pp=setup(root);out=root/'balance';producer=[sys.executable,'scripts/assess_full_screened_balance.py','--plan',str(pp)];reader=[sys.executable,'scripts/readback_full_screened_balance.py','--plan',str(pp),'--output']
        invoke(producer);invoke(reader+[str(root/'passed.json')]);receipt=json.loads((out/'receipt.json').read_text());assert receipt['coverage_rows']==7776 and receipt['balance_rows']==62208
        files=['coverage.tsv','balance.tsv','retained_control_reuse.tsv.gz'];read=lambda name:list(csv.DictReader(gzip.open(out/name,'rt') if name.endswith('.gz') else (out/name).open(),delimiter='\t'));coverage,balance,reuse=[read(name) for name in files];saved={name:(out/name).read_bytes() for name in files};rp=(out/'receipt.json').read_bytes()
        key=('profile',POLICIES[0],'S01','full','n30_c50');cell=next(r for r in coverage if tuple(r[k] for k in KEY)==key);assert cell['retained_matches']=='3' and cell['distinct_control_nodes']=='2' and cell['distinct_control_physical_pairs']=='1' and float(cell['node_weight_kish_concentration'])>2.6 and float(cell['physical_pair_weight_kish_concentration'])==3
        log=next(r for r in balance if tuple(r[k] for k in KEY)==key and r['feature']=='log_positive_sequence_distance');assert log['pairs']=='2' and log['retained_feature_excluded_pairs']=='1'
        constant=next(r for r in balance if tuple(r[k] for k in KEY)==key and r['feature']=='mean_log_length');assert constant['smd_status']=='zero_pooled_variance' and constant['standardized_mean_difference']==''
        assert any(r['scenario_id']=='S54' and r['smd_status']=='insufficient_pairs' for r in balance) and any(r['smd_status']=='no_pairs' for r in balance)
        labels=['missing_coverage','duplicate_coverage','changed_retained_count','changed_taxon_representation','changed_control_node_count','changed_physical_pair_count','changed_weight_sum','changed_weight_concentration','missing_reuse','duplicate_reuse','changed_node_reuse_weight','changed_physical_reuse_weight','changed_log_exclusions','changed_zero_variance_status','changed_original_baseline','changed_matched_retention_shift','changed_eligible_baseline','changed_quantile','changed_summary','changed_scenario54'];rejected=[]
        for label in labels:
            cr,br,rr=copy.deepcopy(coverage),copy.deepcopy(balance),copy.deepcopy(reuse);r=copy.deepcopy(receipt)
            if label=='missing_coverage':cr.pop()
            elif label=='duplicate_coverage':cr.append(copy.deepcopy(cr[0]))
            elif label=='changed_retained_count':cr[0]['retained_matches']='999'
            elif label=='changed_taxon_representation':cr[0]['retained_target_taxa']='999'
            elif label=='changed_control_node_count':cr[0]['distinct_control_nodes']='999'
            elif label=='changed_physical_pair_count':cr[0]['distinct_control_physical_pairs']='999'
            elif label=='changed_weight_sum':cr[0]['reciprocal_node_weight_sum']='999'
            elif label=='changed_weight_concentration':cr[0]['node_weight_kish_concentration']='999'
            elif label=='missing_reuse':rr.pop()
            elif label=='duplicate_reuse':rr.append(copy.deepcopy(rr[0]))
            elif label=='changed_node_reuse_weight':rr[0]['reciprocal_node_reuse_weight']='999'
            elif label=='changed_physical_reuse_weight':rr[0]['reciprocal_physical_pair_reuse_weight']='999'
            elif label=='changed_log_exclusions':next(v for v in br if tuple(v[k] for k in KEY)==key and v['feature']=='log_positive_sequence_distance')['retained_feature_excluded_pairs']='0'
            elif label=='changed_zero_variance_status':next(v for v in br if tuple(v[k] for k in KEY)==key and v['feature']=='mean_log_length')['smd_status']='estimable'
            elif label=='changed_original_baseline':br[0]['baseline_target_mean']='999'
            elif label=='changed_matched_retention_shift':br[0]['metadata_matched_mean_shift']='999'
            elif label=='changed_eligible_baseline':br[0]['target_eligible_all_baseline_targets']='999'
            elif label=='changed_quantile':br[0]['p95_absolute_difference']='999'
            elif label=='changed_summary':r['retained_selection_screen_cells']+=1
            else:next(v for v in rr if v['scenario_id']=='S54')['scenario_id']='S53'
            table(out/files[0],cr);table(out/files[1],br);zipped_table(out/files[2],rr);r['artifacts']={name:sha(out/name) for name in files};write(out/'receipt.json',r);invoke(reader+[str(root/(label+'.json'))],False);rejected.append(label)
        for name,raw in saved.items():(out/name).write_bytes(raw)
        (out/'receipt.json').write_bytes(rp)
        return dict(status='passed_full_screened_balance_software_cases',scenarios=54,strata=432,coverage_rows=7776,balance_rows=62208,selected_records=8,three_baselines_checked=True,zero_log_distance_exclusions_and_constant_variance_checked=True,node_and_physical_pair_reuse_distinguished=True,all_masks_screens_and_empty_strata_checked=True,rejected_rehashed_exports=rejected,synthetic_prior_matching_physics_and_journal_contracts=True,production_physics_or_journals_tested=False,scope=__doc__)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run();result['script_sha256']=sha(__file__)
    with args.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
