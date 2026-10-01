#!/usr/bin/env python3
"""Complete 54-scenario software fixture, preserved exclusions and rehashed false exports.

All prior physical coverage/native/matching/journal proof records are synthetic
contracts. Full source I/O, node/physical-pair role projection, mask intersections,
frozen selection/unmatched membership, export fidelity and independent SQL
reconstruction run unchanged. This is not physical validation, real matching
proof, production closure or a biological pilot.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

SCREENS=[dict(id=f'n{n}_c{c}',minimum_aligned_residues=n,minimum_original_coverage=c/100) for n in [30,50] for c in [50,70,90]]
POLICIES=['alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore']


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')


def table(path,rows):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def zipped_table(path,rows):
    with gzip.open(path,'wt') as handle:
        writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def invoke(command,success=True):
    result=subprocess.run(command,capture_output=True,text=True)
    if success and result.returncode:raise RuntimeError(result.stderr+result.stdout)
    if not success:assert result.returncode!=0,'False export accepted'


def node(ident,guide,a,b,target=False,reverse=False):
    ends=[(a,6),(b,10 if b in ['B','G','J'] else 6)]
    if a==b:ends=[(a,6),(a,6)]
    if reverse:ends.reverse()
    lengths={a:100,b:100 if a==b else 101}
    value=dict(node_id=ident,guide=guide,family='OG-fixture',gene_a=ident+'-gene-a',gene_b=ident+'-gene-b',sequence_distance=0.25,pair_key=hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest(),same_model=int(ends[0]==ends[1]))
    for side,end in zip(['a','b'],ends):value.update({f'model_id_{side}':end[0],f'version_{side}':end[1],f'length_{side}':lengths[end[0]]})
    if target:value.update(taxon_id='taxon-'+guide,gene_node='node-'+ident)
    else:value.update(taxon_a='taxon-'+guide,taxon_b='taxon-other',both_guides=1,both_unreported_parents=1)
    return value


def quality_rows(nodes,bits):
    records=[]
    for n,(full,masked) in zip(nodes,bits):
        if n['same_model']:continue
        sides=sorted([(n['model_id_'+s],n['version_'+s],n['length_'+s]) for s in ['a','b']])
        for mask,value in [('full',full),('plddt70',masked)]:
            row=dict(pair_key=n['pair_key'],mask=mask,model_a=sides[0][0],version_a=sides[0][1],model_b=sides[1][0],version_b=sides[1][1],length_a=sides[0][2],length_b=sides[1][2])
            for j,spec in enumerate(SCREENS):row.update({spec['id']+'_pass':int(bool(value&(1<<j))),spec['id']+'_exclusions':'' if value&(1<<j) else 'synthetic_retained_exclusion'})
            records.append(row)
    return records


def setup(root):
    target,background,selection,graph=[root/name for name in ['target','background','selection','graph']]
    for p in [target,background,selection,graph]:p.mkdir()
    targets=[node('T1','profile','A','B',True),node('T2','profile','A','A',True),node('T3','mafft','A','C',True,True),node('T4','mafft','D','E',True)]
    backgrounds=[node('B1','profile','F','G'),node('B2','profile','F','H'),node('B3','mafft','I','I'),node('B4','mafft','J','K')]
    tq=quality_rows(targets,[(63,9),(0,0),(63,0),(0,0)]);bq=quality_rows(backgrounds,[(63,9),(0,0),(0,0),(63,63)])
    table(target/'pair_mask_coverage.tsv',tq);table(background/'pair_mask_coverage.tsv',bq)
    pass_counts=lambda rows:{mask+':'+s['id']:sum(r[s['id']+'_pass'] for r in rows if r['mask']==mask) for mask in ['full','plddt70'] for s in SCREENS}
    both_counts=lambda rows:{s['id']:sum(all(next(r for r in rows if r['pair_key']==key and r['mask']==mask)[s['id']+'_pass'] for mask in ['full','plddt70']) for key in {r['pair_key'] for r in rows}) for s in SCREENS}
    tp,bp=root/'target-plan.json',root/'background-plan.json';write(tp,dict(output=str(target),screens=SCREENS));write(bp,dict(output=str(background)))
    tr=target/'receipt.json';t=dict(status='complete_full_expanded_duplication_coverage_pending_independent_readback',plan_sha256=sha(tp),pairs=3,pair_mask_rows=6,screens=SCREENS,pair_pass_counts=pass_counts(tq),pair_both_masks_pass_counts=both_counts(tq),artifacts={'pair_mask_coverage.tsv':sha(target/'pair_mask_coverage.tsv')});write(tr,t)
    tc=root/'target-completed.json';write(tc,dict(status='complete_verified_full_expanded_pair_event_and_taxon_coverage_screens',services=[dict(synthetic_fixture_only=True)]*2,scientific_eligibility=False,summary={k:t[k] for k in ['pairs','pair_mask_rows','pair_pass_counts','pair_both_masks_pass_counts']},source_hashes={str(p):sha(p) for p in [tp,tr,target/'pair_mask_coverage.tsv']}))
    summary=dict(source_inventory_models=8,pair_models=5,pairs=3,pair_mask_rows=6,directed_source_states=12,screen_decisions=36,screens=SCREENS,pass_counts=pass_counts(bq),both_masks_pass_counts=both_counts(bq),exclusion_counts={},target_pairs=3,target_pass_counts=t['pair_pass_counts'],target_both_masks_pass_counts=t['pair_both_masks_pass_counts'])
    br,ba=background/'receipt.json',background/'readback.json';write(br,dict(status='complete_full_background_coverage_pending_independent_readback',plan_sha256=sha(bp),**summary,scientific_eligibility=False,artifacts={'pair_mask_coverage.tsv':sha(background/'pair_mask_coverage.tsv')}));write(ba,dict(status='passed_full_background_coverage_sql_decimal_readback',plan_sha256=sha(bp),producer_receipt_sha256=sha(br),**summary,scientific_eligibility=False))
    archive,bc=background/'archive.json',root/'background-completed.json';bindings={str(p):sha(p) for p in [bp,br,ba,background/'pair_mask_coverage.tsv']};write(archive,dict(status='complete_verified_full_background_coverage_archive',summary=summary,services=[dict(synthetic_fixture_only=True)]*2,source_hashes=bindings))
    write(bc,dict(status='complete_verified_full_background_original_length_coverage',**summary,source_plan=str(bp),source_plan_sha256=sha(bp),producer_receipt=str(br),producer_receipt_sha256=sha(br),independent_readback=str(ba),independent_readback_sha256=sha(ba),scientific_eligibility=False,exact_process_journals_checked=2,bound_source_hashes=len(bindings),full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive)))
    for name,records in [('target_nodes.jsonl',targets),('background_nodes.jsonl',backgrounds)]:
        (graph/name).write_text(''.join(json.dumps(n)+'\n' for n in records))
    write(graph/'receipt.json',dict(status='complete_background_match_graph_pending_readback',target_nodes=4,background_nodes=4,target_policy_dispositions=16,edges=6,artifacts={name:sha(graph/name) for name in ['target_nodes.jsonl','background_nodes.jsonl']}))
    covariates=root/'covariates';covariates.mkdir();gp=graph/'receipt.json';cp=covariates/'receipt.json';gap=root/'graph-readback.json';cap=root/'covariate-readback.json';gc=root/'graph-completed.json'
    write(gap,dict(status='passed_full_background_match_graph_readback',producer_receipt_sha256=sha(gp),target_nodes=4,background_nodes=4,target_policy_dispositions=16,edges=6,support_rows_checked=48))
    write(cp,dict(status='complete_background_match_covariates_pending_readback',graph_receipt_sha256=sha(gp),graph_readback_sha256=sha(gap),edges=6))
    write(cap,dict(status='passed_full_background_match_covariate_readback',producer_receipt_sha256=sha(cp),edges=6))
    write(gc,dict(status='complete_verified_expanded_background_graph_covariates',services=[dict(synthetic_fixture_only=True)]*4,scientific_eligibility=False,summary=dict(duplicate_target_links=4,background_nodes=4,target_policy_dispositions=16,eligible_edges=6,architecture_support_rows_checked=48),source_hashes={str(p):sha(p) for p in [gp,cp,gap,cap,graph/'target_nodes.jsonl',graph/'background_nodes.jsonl']}))
    scenarios=[dict(scenario_id=f'S{i:02d}',fixture_only=True) for i in range(1,55)];write(selection/'scenarios.json',scenarios)
    choices={('T1',POLICIES[0]):[('S01','B1'),('S54','B1')],('T2',POLICIES[0]):[('S01','B2')],('T3',POLICIES[0]):[('S01','B3')],('T4',POLICIES[0]):[('S01','B4')],('T1',POLICIES[1]):[('S01','B2')]}
    statuses=[];chosen=[];counts=Counter()
    for n in targets:
        for p in POLICIES:
            used=choices.get((n['node_id'],p),[]);matched={sid for sid,bid in used}
            statuses.append(dict(target_id=n['node_id'],policy=p,architecture_status='fixture',candidate_edges=len(used),shared_identity_edges=0,matched_scenarios=','.join(s['scenario_id'] for s in scenarios if s['scenario_id'] in matched),unmatched_scenarios=','.join(s['scenario_id'] for s in scenarios if s['scenario_id'] not in matched)))
            for s in scenarios:
                prefix='|'.join([n['guide'],p,s['scenario_id']]);counts[prefix+'|targets']+=1;counts[prefix+('|matched' if s['scenario_id'] in matched else '|unmatched')]+=1
            for sid,bid in used:chosen.append(dict(target_id=n['node_id'],policy=p,scenario_id=sid,background_id=bid,endpoint_order='1',score='0.25',eligible_candidates='2',equal_score_candidates='2',score_gap_to_second='0.0'))
    table(selection/'target_policy_selection_status.tsv',statuses);zipped_table(selection/'selections.tsv.gz',chosen)
    mr=selection/'receipt.json';matching_summary=dict(target_policy_records=16,scenarios=54,scenario_decisions=864,selected_records=6,unmatched_decisions=858,source_edges=6)
    write(mr,dict(status='complete_metadata_background_control_selection_pending_readback',source_covariate_receipt_sha256=sha(cp),**matching_summary,counts=dict(counts),artifacts={name:sha(selection/name) for name in ['scenarios.json','target_policy_selection_status.tsv','selections.tsv.gz']}))
    mc=root/'matching-completed.json';paths=[mr,selection/'scenarios.json',selection/'target_policy_selection_status.tsv',selection/'selections.tsv.gz'];write(mc,dict(status='complete_verified_expanded_background_fixed_matching',services=[dict(synthetic_fixture_only=True)]*2,scientific_eligibility=False,summary=matching_summary,source_hashes={str(p):sha(p) for p in paths}))
    census_plan=root/'census-plan.json';census=root/'census.json';census_completion=root/'census-completed.json'
    write(census_plan,dict(synthetic_fixture_only=True));census_summary=dict(target_nodes=4,background_nodes=4,target_physical_pairs=3,background_physical_pairs=3,source_inventory_models=8,used_target_pairs=3,used_background_pairs=3,confidence_descriptor_mismatch_endpoint_occurrences={},maximum_confidence_descriptor_absolute_differences={})
    write(census,dict(status='passed_complete_fixed_matching_node_physical_identity_census',plan_sha256=sha(census_plan),**census_summary,source_hashes={str(p):sha(p) for p in [graph/'receipt.json',graph/'target_nodes.jsonl',graph/'background_nodes.jsonl']}))
    write(census_completion,dict(status='complete_verified_full_fixed_matching_node_physical_identity_census',services=[dict(synthetic_fixture_only=True)],scientific_eligibility=False,summary=census_summary,source_hashes={str(p):sha(p) for p in [census_plan,census,graph/'receipt.json',graph/'target_nodes.jsonl',graph/'background_nodes.jsonl']}))
    pp=root/'plan.json';write(pp,dict(background_completion=str(bc),background_plan=str(bp),target_completion=str(tc),target_plan=str(tp),matching_completion=str(mc),graph_completion=str(gc),covariates=str(covariates),graph_readback=str(gap),covariate_readback=str(cap),selection=str(selection),graph=str(graph),screens=SCREENS,guides=['profile','mafft'],policies=POLICIES,output=str(root/'integrated'),expected=dict(background_pairs=3,target_pairs=3,target_nodes=4,background_nodes=4,target_policy_records=16,scenarios=54,selected_records=6,unmatched_decisions=858,source_edges=6),pins={},resources=dict(minimum_free_disk_gib=0),scope=__doc__))
    config=json.loads(pp.read_text());config.update(input_census=str(census),input_census_completion=str(census_completion),input_census_plan=str(census_plan));config['expected']['source_inventory_models']=8;write(pp,config)
    return pp


def run():
    with tempfile.TemporaryDirectory(prefix='full-matched-coverage-fixture-') as temp:
        root=Path(temp);pp=setup(root);out=root/'integrated';producer=[sys.executable,'scripts/integrate_full_matched_coverage.py','--plan',str(pp)];reader=[sys.executable,'scripts/readback_full_matched_coverage.py','--plan',str(pp),'--output']
        invoke(producer);invoke(reader+[str(root/'passed.json')])
        receipt=json.loads((out/'receipt.json').read_text());assert receipt['selected_records']==6 and receipt['unmatched_decisions']==858 and receipt['attrition_rows']==7776 and receipt['full_scenario_screen_cells']==15552
        files=['selected_pair_coverage.tsv.gz','target_policy_coverage_status.tsv.gz','matched_attrition_counts.tsv'];saved={name:(out/name).read_bytes() for name in files};rp=(out/'receipt.json').read_bytes()
        read=lambda path:list(csv.DictReader(gzip.open(path,'rt') if str(path).endswith('.gz') else Path(path).open(),delimiter='\t'))
        selections=read(out/files[0]);status=read(out/files[1]);attrition=read(out/files[2]);lookup={tuple(row[k] for k in ['guide','policy','scenario_id','mask','screen']):row for row in attrition}
        assert lookup['profile',POLICIES[0],'S01','full','n30_c50']['joint_pass_matched_records']=='1'
        assert lookup['mafft',POLICIES[0],'S01','full','n30_c50']['target_only_pass_matched_records']=='1' and lookup['mafft',POLICIES[0],'S01','full','n30_c50']['control_only_pass_matched_records']=='1'
        assert next(r for r in selections if r['scenario_id']=='S54')['joint_plddt70_pass_bits']=='9'
        labels=['missing_selection','duplicate_selection','changed_fixed_control','changed_endpoint_order','changed_match_score','changed_target_bits','changed_control_bits','changed_joint_bits','changed_pair_identity','promoted_same_model','missing_policy','changed_unmatched_scenarios','changed_high_scenario','changed_summary','changed_attrition','missing_empty_stratum','duplicate_stratum'];rejected=[]
        for label in labels:
            sr,ur,ar=copy.deepcopy(selections),copy.deepcopy(status),copy.deepcopy(attrition);r=copy.deepcopy(receipt)
            if label=='missing_selection':sr.pop()
            elif label=='duplicate_selection':sr.append(copy.deepcopy(sr[0]))
            elif label=='changed_fixed_control':sr[0]['background_id']='B2'
            elif label=='changed_endpoint_order':sr[0]['endpoint_order']='0'
            elif label=='changed_match_score':sr[0]['score']='999'
            elif label=='changed_target_bits':sr[0]['target_full_pass_bits']='0'
            elif label=='changed_control_bits':sr[0]['control_plddt70_pass_bits']='63'
            elif label=='changed_joint_bits':sr[0]['joint_both_pass_bits']='63'
            elif label=='changed_pair_identity':sr[0]['target_pair_key']='0'*64
            elif label=='promoted_same_model':next(row for row in sr if row['target_id']=='T2')['target_full_pass_bits']='63'
            elif label=='missing_policy':ur.pop()
            elif label=='changed_unmatched_scenarios':ur[0]['unmatched_scenarios']=''
            elif label=='changed_high_scenario':next(row for row in sr if row['scenario_id']=='S54')['scenario_id']='S53'
            elif label=='changed_summary':r['unmatched_decisions']-=1
            elif label=='changed_attrition':ar[0]['joint_pass_matched_records']='999'
            elif label=='missing_empty_stratum':ar.pop(next(i for i,row in enumerate(ar) if row['matched_records']=='0'))
            else:ar.append(copy.deepcopy(ar[0]))
            zipped_table(out/files[0],sr);zipped_table(out/files[1],ur);table(out/files[2],ar);r['artifacts']={name:sha(out/name) for name in files};write(out/'receipt.json',r)
            invoke(reader+[str(root/(label+'.json'))],success=False);rejected.append(label)
        for name,blob in saved.items():(out/name).write_bytes(blob)
        (out/'receipt.json').write_bytes(rp)
        return dict(status='passed_full_matched_coverage_software_cases',guides=2,policies=4,scenarios=54,target_policy_records=16,selected_records=6,unmatched_decisions=858,full_scenario_screen_cells=15552,attrition_rows=7776,reversed_gene_model_roles_and_versions6_10_checked=True,same_model_exclusions_retained=True,all_four_attrition_categories_checked=True,scenario54_and_empty_strata_checked=True,rejected_rehashed_exports=rejected,synthetic_prior_physical_matching_and_journal_contracts=True,production_native_matching_or_journal_correctness_tested=False,scope=__doc__)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args();result=run();result['script_sha256']=sha(__file__)
    if args.output:
        with args.output.open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
