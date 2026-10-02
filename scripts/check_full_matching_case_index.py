#!/usr/bin/env python3
"""Exercise full-design indexing, distinct genes/shared models and false exports.

Prior measurement, matching and journal proofs are synthetic contracts. These
checks verify software identity/membership/census behavior, not physical
measurements, full scientific source acceptance or a biological pilot.
"""
import argparse
from collections import Counter
import copy
import csv
import gzip
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from check_full_matched_coverage_cases import node,write,table,zipped_table,SCREENS,POLICIES
from full_matching_case_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def invoke(cmd, success=True):
    result=subprocess.run(cmd,text=True,capture_output=True)
    if success and result.returncode:raise RuntimeError(result.stderr+result.stdout)
    if not success:assert result.returncode!=0,'Invalid case index accepted'


def setup(root):
    graph,matching,selection=[root/n for n in ['graph','matching','selection']]
    for p in [graph,matching,selection]:p.mkdir()
    targets=[node('T1','profile','A','B',True),node('T2','profile','A','B',True,True),node('T3','profile','A','A',True),node('T4','mafft','C','D',True),node('T5','mafft','E','F',True)]
    backgrounds=[node('B1','profile','G','H'),node('B2','profile','G','H',reverse=True),node('B3','mafft','I','I'),node('B4','mafft','J','K')]
    ts={n['node_id']:n for n in targets};bs={n['node_id']:n for n in backgrounds}
    for kind,ns in [('target',targets),('background',backgrounds)]: (graph/(kind+'_nodes.jsonl')).write_text(''.join(json.dumps(n)+'\n' for n in ns))
    tf={'T1':[63,9,9],'T2':[63,9,9],'T3':[0,0,0],'T4':[9,0,0],'T5':[0,0,0]};bf={'B1':[63,9,9],'B2':[63,9,9],'B3':[0,0,0],'B4':[0,0,0]}
    scenarios=[dict(scenario_id=f'S{i:02d}',synthetic_fixture_only=True) for i in range(1,55)];write(selection/'scenarios.json',scenarios)
    choices=[('T1',POLICIES[0],'S01','B1',0),('T1',POLICIES[0],'S54','B1',1),('T1',POLICIES[1],'S01','B2',1),('T2',POLICIES[0],'S01','B1',1),('T3',POLICIES[0],'S01','B2',0),('T4',POLICIES[0],'S01','B3',1),('T4',POLICIES[1],'S54','B4',0)]
    selected=[];statuses=[]
    for t,p in itertools.product(targets,POLICIES):
        matched=[]
        for tid,policy,sid,bid,order in choices:
            if (tid,policy)!=(t['node_id'],p):continue
            b=bs[bid];matched.append(sid)
            row=dict(target_id=tid,policy=policy,scenario_id=sid,background_id=bid,endpoint_order=order,score=1.5,eligible_candidates=2,equal_score_candidates=1,score_gap_to_second='',guide=t['guide'],family=t['family'],focal_taxon=t['taxon_id'],gene_node=t['gene_node'],target_pair_key=t['pair_key'],background_pair_key=b['pair_key'],target_comparison_disposition='identical_model_no_alignment' if t['same_model'] else 'distinct_model_pair',background_comparison_disposition='identical_model_no_alignment' if b['same_model'] else 'distinct_model_pair',target_sequence_distance=t['sequence_distance'],background_sequence_distance=b['sequence_distance'])
            for i,m in enumerate(['full','plddt70','both']):row.update({f'target_{m}_pass_bits':tf[tid][i],f'control_{m}_pass_bits':bf[bid][i],f'joint_{m}_pass_bits':tf[tid][i]&bf[bid][i]})
            selected.append(row)
        statuses.append(dict(target_id=t['node_id'],policy=p,architecture_status='synthetic_only',candidate_edges=2,shared_identity_edges=2,matched_scenarios=','.join(matched),unmatched_scenarios=','.join(s['scenario_id'] for s in scenarios if s['scenario_id'] not in matched),guide=t['guide'],family=t['family'],focal_taxon=t['taxon_id'],gene_node=t['gene_node'],target_pair_key=t['pair_key'],target_comparison_disposition='identical_model_no_alignment' if t['same_model'] else 'distinct_model_pair',**{f'target_{m}_pass_bits':tf[t['node_id']][i] for i,m in enumerate(['full','plddt70','both'])}))
    zipped_table(matching/'selected_pair_coverage.tsv.gz',selected);zipped_table(matching/'target_policy_coverage_status.tsv.gz',statuses)
    attrs=[]
    for g,p,s,m,spec in itertools.product(['profile','mafft'],POLICIES,scenarios,['full','plddt70','both'],SCREENS):
        i=['full','plddt70','both'].index(m);bit=SCREENS.index(spec);ns=[n for n in targets if n['guide']==g];chosen=[r for r in selected if (r['guide'],r['policy'],r['scenario_id'])==(g,p,s['scenario_id'])]
        cats=Counter((int(bool(tf[r['target_id']][i]&1<<bit)),int(bool(bf[r['background_id']][i]&1<<bit))) for r in chosen);tp=cats[1,0]+cats[1,1];ap=sum(bool(tf[n['node_id']][i]&1<<bit) for n in ns)
        attrs.append(dict(guide=g,policy=p,scenario_id=s['scenario_id'],mask=m,screen=spec['id'],all_target_records=len(ns),matched_records=len(chosen),unmatched_records=len(ns)-len(chosen),target_pass_all_records=ap,target_pass_matched_records=tp,target_pass_unmatched_records=ap-tp,control_pass_matched_records=cats[0,1]+cats[1,1],joint_pass_matched_records=cats[1,1],target_only_pass_matched_records=cats[1,0],control_only_pass_matched_records=cats[0,1],neither_pass_matched_records=cats[0,0]))
    table(matching/'matched_attrition_counts.tsv',attrs)
    expected=dict(target_nodes=5,background_nodes=4,selected_records=7,scenarios=54,target_policy_records=20,unmatched_decisions=1073)
    mp=root/'matching-plan.json';write(mp,dict(output=str(matching),graph=str(graph),selection=str(selection),screens=SCREENS,expected=expected))
    summary=dict(**expected,scenario_decisions=1080,selection_screen_cells=126,full_scenario_screen_cells=19440,attrition_rows=7776,screens=SCREENS,masks=['full','plddt70','both'],guides=['profile','mafft'],policies=POLICIES,node_dispositions={})
    rp,ap=matching/'receipt.json',matching/'readback.json';write(rp,dict(status='complete_full_matched_coverage_pending_independent_readback',plan_sha256=sha(mp),**summary,artifacts={n:sha(matching/n) for n in ['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv']}));write(ap,dict(status='passed_full_matched_coverage_sql_readback',plan_sha256=sha(mp),producer_receipt_sha256=sha(rp),**summary))
    paths=[mp,rp,ap,graph/'target_nodes.jsonl',graph/'background_nodes.jsonl',selection/'scenarios.json']+[matching/n for n in ['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv']];bindings={str(p):sha(p) for p in paths}
    archive=matching/'archive.json';write(archive,dict(status='complete_verified_full_matched_coverage_archive',summary=summary,services=[dict(synthetic_fixture_only=True)]*2,source_hashes=bindings))
    cp=root/'matching-completed.json';write(cp,dict(status='complete_verified_full_fixed_matched_coverage_attrition',**summary,source_plan=str(mp),source_plan_sha256=sha(mp),producer_receipt=str(rp),producer_receipt_sha256=sha(rp),independent_readback=str(ap),independent_readback_sha256=sha(ap),full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),bound_source_hashes=len(bindings),exact_process_journals_checked=2,scientific_eligibility=False))
    pp=root/'plan.json';write(pp,dict(matching_completion=str(cp),matching_plan=str(mp),screens=SCREENS,guides=['profile','mafft'],policies=POLICIES,expected=expected,output=str(root/'cases'),pins={},resources=dict(minimum_free_disk_gib=0),scope=__doc__));return pp


def run():
    with tempfile.TemporaryDirectory(prefix='full-matching-case-fixture-') as temp:
        root=Path(temp);pp=setup(root);out=root/'cases';producer=[sys.executable,'scripts/index_full_matching_cases.py','--plan',str(pp)];reader=[sys.executable,'scripts/readback_full_matching_cases.py','--plan',str(pp),'--output']
        invoke(producer);invoke(reader+[str(root/'passed.json')]);receipt=json.loads((out/'receipt.json').read_text())
        assert receipt['logical_cases']==6 and receipt['physical_cases']==4 and receipt['selected_records']==7 and receipt['unmatched_decisions']==1073
        read=lambda name:list(csv.DictReader(gzip.open(out/name,'rt'),delimiter='\t'))
        cases=read('case_index.tsv.gz');links=read('selection_case_links.tsv.gz')
        c=next(r for r in cases if (r['target_id'],r['background_id'])==('T1','B1'))
        assert c['endpoint_order_bits']=='3' and c['scenario_bits']==str(1+2**53) and c['selection_records']=='2'
        assert len({r['case_id'] for r in cases if r['physical_case_id']==c['physical_case_id']})==3
        assert any(r['target_same_model']=='1' for r in cases) and any(r['background_same_model']=='1' for r in cases)
        invoke(producer,False)
        files=['case_index.tsv.gz','selection_case_links.tsv.gz'];saved={n:(out/n).read_bytes() for n in files};original=(out/'receipt.json').read_bytes()
        labels=['missing_case','duplicate_case','collapsed_logical_id','wrong_physical_id','wrong_gene_node','wrong_gene_identity','wrong_family','wrong_model_pair','wrong_case_reuse','wrong_endpoint_bits','wrong_policy_bits','missing_scenario54','promoted_same_model','wrong_joint_mask','missing_selection','duplicate_selection','wrong_selection_ordinal','changed_endpoint_mapping','changed_fixed_score','changed_fixed_control','wrong_selection_link','changed_summary','changed_unmatched_locator']
        for label in labels:
            cr,lr=copy.deepcopy(cases),copy.deepcopy(links);r=copy.deepcopy(receipt)
            if label=='missing_case':cr.pop()
            elif label=='duplicate_case':cr.append(copy.deepcopy(cr[0]))
            elif label=='collapsed_logical_id':cr[0]['case_id']=cr[0]['physical_case_id']
            elif label=='wrong_physical_id':cr[0]['physical_case_id']='0'*64
            elif label=='wrong_gene_node':cr[0]['gene_node']='other'
            elif label=='wrong_gene_identity':cr[0]['background_gene_a']='other'
            elif label=='wrong_family':cr[0]['background_family']='other'
            elif label=='wrong_model_pair':cr[0]['target_pair_key']='0'*64
            elif label=='wrong_case_reuse':cr[0]['selection_records']='1'
            elif label=='wrong_endpoint_bits':cr[0]['endpoint_order_bits']='1'
            elif label=='wrong_policy_bits':cr[0]['policy_bits']='15'
            elif label=='missing_scenario54':cr[0]['scenario_bits']='1'
            elif label=='promoted_same_model':next(c for c in cr if c['target_same_model']=='1')['target_full_pass_bits']='63'
            elif label=='wrong_joint_mask':cr[0]['joint_both_pass_bits']='63'
            elif label=='missing_selection':lr.pop()
            elif label=='duplicate_selection':lr.append(copy.deepcopy(lr[0]))
            elif label=='wrong_selection_ordinal':lr[0]['source_row_ordinal']='2'
            elif label=='changed_endpoint_mapping':lr[0]['endpoint_order']='1'
            elif label=='changed_fixed_score':lr[0]['score']='999'
            elif label=='changed_fixed_control':lr[0]['background_id']='B2'
            elif label=='wrong_selection_link':lr[0]['case_id']='0'*64
            elif label=='changed_summary':r['physical_cases']+=1
            else:r['original_unmatched_status_sha256']='0'*64
            zipped_table(out/files[0],cr);zipped_table(out/files[1],lr);r['artifacts'].update({n:sha(out/n) for n in files});write(out/'receipt.json',r)
            invoke(reader+[str(root/(label+'.json'))],False)
        for n,b in saved.items():(out/n).write_bytes(b)
        (out/'receipt.json').write_bytes(original)
        # Simulate an interrupted, unaccepted artifact publication. The original
        # plan is required, a full replay replaces scratch, and counts survive.
        (out/'receipt.json').unlink();(out/'selection_case_links.tsv.gz.tmp').write_text('incomplete')
        invoke(producer);invoke(reader+[str(root/'recovered.json')]);recovered=json.loads((out/'receipt.json').read_text())
        assert all(receipt[k]==recovered[k] for k in SUMMARY_FIELDS) and cases==read(files[0]) and links==read(files[1])
        return dict(status='passed_full_matching_case_index_software_contracts',selected_records=7,logical_cases=6,physical_cases=4,scenarios=54,unmatched_decisions=1073,distinct_genes_sharing_models_preserved=True,both_endpoint_mappings_and_scenario54_checked=True,same_model_quality_exclusions_preserved=True,completed_restart_refused=True,interrupted_full_replay_passed=True,rejected_rehashed_exports=labels,synthetic_prior_measurement_matching_and_journal_contracts=True,production_physics_or_journals_tested=False,scope=__doc__)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run();r['source_hashes']={n:sha(n) for n in ['scripts/full_matching_case_sources.py','scripts/index_full_matching_cases.py','scripts/readback_full_matching_cases.py',__file__]}
    with a.output.open('x') as handle:handle.write(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
