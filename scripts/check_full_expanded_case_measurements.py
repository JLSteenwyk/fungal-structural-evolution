#!/usr/bin/env python3
"""Complete case/mask/order-contrast join and independent Decimal contract checks.

Prior native fits, fixed matching, indexing and journals are synthetic source
contracts. No biological subset or real scientific acceptance is claimed.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from check_full_expanded_measurement_catalog import setup as catalog_setup,invoke
from check_full_matched_coverage_cases import write,zipped_table
from full_matching_case_sources import CASE_FIELDS,identity
from full_expanded_case_measurement_sources import SUMMARY_FIELDS
from run_ortholog_pair_guide_comparison import sha


def close(root,plan,completed_status,receipt):
    rp=root/'receipt.json';readback=root/'synthetic_prior_readback.json'
    write(readback,dict(status='synthetic_prior_readback_contract',plan_sha256=sha(plan),producer_receipt_sha256=sha(rp),scientific_eligibility=False))
    bindings=dict(receipt.get('source_hashes',{}));bindings.update({str(p):sha(p) for p in [plan,rp,readback]})
    for n,d in receipt['artifacts'].items():bindings[str(root/n)]=d
    summary={k:receipt[k] for k in ['selected_records','logical_cases','physical_cases','unmatched_decisions','directed_states'] if k in receipt}
    archive=root/'synthetic_completion_archive.json';write(archive,dict(status=completed_status+'_archive',services=[dict(synthetic_fixture_only=True)]*2,source_hashes=bindings,summary=summary))
    c=root/'synthetic_completed.json';write(c,dict(status=completed_status,**summary,scientific_eligibility=False,exact_process_journals_checked=2,bound_source_hashes=len(bindings),full_hash_archive=str(archive),full_hash_archive_sha256=sha(archive),producer_receipt=str(rp),producer_receipt_sha256=sha(rp)));return c


def setup(root):
    cp=catalog_setup(root);invoke([sys.executable,'scripts/export_full_expanded_measurement_catalog.py','--plan',str(cp)])
    config=json.loads(cp.read_text());catalog=Path(config['output']);cr=json.loads((catalog/'receipt.json').read_text());catalog_closed=close(catalog,cp,'complete_verified_full_expanded_directed_measurement_catalog',cr)
    states=list(csv.DictReader(gzip.open(catalog/'directed_measurements.tsv.gz','rt'),delimiter='\t'))
    pairs=sorted({r['pair_key'] for r in states});same=identity('synthetic-same-model-pair','same','same')
    qs={(r['role'],r['pair_key'],r['mask']):r for r in states}
    layouts=[(pairs[0],pairs[0]),(pairs[0],pairs[1]),(same,pairs[1]),(pairs[1],same),(pairs[0],pairs[0]),(pairs[1],pairs[1])]
    cases=[]
    for i,(tp,bp) in enumerate(layouts):
        tid,bid=f'T{i}',f'B{i}';cid=identity('fixed-matched-logical-case-v1',tid,bid)
        row={k:'' for k in CASE_FIELDS};row.update(case_id=cid,physical_case_id=identity('fixed-matched-physical-case-v1',tp,bp),target_id=tid,background_id=bid,guide='profile',target_family='OG-target',background_family='OG-control',focal_taxon='taxon',gene_node='gene-node-'+tid,target_gene_a=tid+'-a',target_gene_b=tid+'-b',background_gene_a=bid+'-a',background_gene_b=bid+'-b',background_taxon_a='taxon-a',background_taxon_b='taxon-b',target_pair_key=tp,background_pair_key=bp,target_same_model=int(tp==same),background_same_model=int(bp==same),target_sequence_distance=0.0 if i==0 else .25,background_sequence_distance=.1,selection_records=2 if i==0 else 1,endpoint_order_bits=3 if i==0 else 1,policy_bits=1,scenario_bits=1)
        for mask in ['full','plddt70']:
            tf=0 if tp==same else int(qs['target',tp,mask]['pair_mask_pass_bits']);bf=0 if bp==same else int(qs['background',bp,mask]['pair_mask_pass_bits']);row.update({f'target_{mask}_pass_bits':tf,f'control_{mask}_pass_bits':bf,f'joint_{mask}_pass_bits':tf&bf})
        for role in ['target','control','joint']:row[role+'_both_pass_bits']=row[role+'_full_pass_bits']&row[role+'_plddt70_pass_bits']
        cases.append(row)
    index=root/'case-index';index.mkdir();zipped_table(index/'case_index.tsv.gz',cases)
    ip=root/'case-index-plan.json';write(ip,dict(output=str(index),matching_completion=config['matching_completion'],matching_plan=config['matching_plan']))
    ir=dict(status='complete_full_matching_case_index_pending_independent_readback',plan_sha256=sha(ip),scientific_eligibility=False,selected_records=7,logical_cases=6,physical_cases=5,unmatched_decisions=1073,source_hashes={},artifacts={'case_index.tsv.gz':sha(index/'case_index.tsv.gz')});write(index/'receipt.json',ir)
    index_closed=close(index,ip,'complete_verified_full_matching_logical_case_index',ir)
    pp=root/'join-plan.json';write(pp,dict(case_index_completion=str(index_closed),case_index_plan=str(ip),catalog_completion=str(catalog_closed),catalog_plan=str(cp),matching_completion=config['matching_completion'],expected=dict(selected_records=7,logical_cases=6,unmatched_decisions=1073,directed_states=16),pins={},resources=dict(minimum_free_disk_gib=0),output=str(root/'joined'),scope=__doc__));return pp


def run():
    with tempfile.TemporaryDirectory(prefix='expanded-case-measurements-fixture-') as temp:
        root=Path(temp);pp=setup(root);out=root/'joined';producer=[sys.executable,'scripts/join_full_expanded_case_measurements.py','--plan',str(pp)];reader=[sys.executable,'scripts/readback_full_expanded_case_measurements.py','--plan',str(pp),'--output']
        invoke(producer);invoke(reader+[str(root/'passed.json')]);rp=out/'receipt.json';r=json.loads(rp.read_text());assert r['case_mask_rows']==12 and r['order_pair_cells']==96
        data=out/'case_mask_contrasts.tsv.gz';rows=list(csv.DictReader(gzip.open(data,'rt'),delimiter='\t'))
        complete=next(x for x in rows if x['rmsd_usable_order_pair_bits']=='15');incomplete=next(x for x in rows if x['rmsd_usable_order_pair_bits']!='15')
        assert complete['rmsd_complete_order_envelope_min']=='-1.0' and complete['rmsd_complete_order_envelope_max']=='1.0' and complete['rmsd_both_orders_mean_delta']=='0.0'
        assert incomplete['rmsd_complete_order_envelope_min']==''
        invoke(producer,False);saved=data.read_bytes();original=rp.read_bytes()
        labels=['missing_row','duplicate_row','changed_row_identity','collapsed_gene_case','changed_mask','changed_sequence_distance','changed_state_reference','changed_native_failure','quarantined_rmsd_promoted','missing_as_zero','changed_order_pair_01','favorable_order_mean','partial_envelope_promoted','changed_envelope_span','changed_pair_bitmap','changed_native_tm_dissimilarity','changed_joint_quality','changed_same_model','changed_summary']
        for label in labels:
            changed=copy.deepcopy(rows);receipt=copy.deepcopy(r)
            if label=='missing_row':changed.pop()
            elif label=='duplicate_row':changed.append(copy.deepcopy(changed[0]))
            elif label=='changed_row_identity':changed[0]['row_identity']='0'*64
            elif label=='collapsed_gene_case':changed[0]['case_id']=changed[0]['physical_case_id']
            elif label=='changed_mask':changed[0]['mask']='both'
            elif label=='changed_sequence_distance':changed[0]['target_sequence_distance']='999'
            elif label=='changed_state_reference':changed[0]['target_order0_state_key']='other'
            elif label=='changed_native_failure':changed[0]['target_order0_native_status']='timeout'
            elif label=='quarantined_rmsd_promoted':next(x for x in changed if x['target_order0_numerical_usable']=='0')['target_order0_rmsd']='1.0'
            elif label=='missing_as_zero':next(x for x in changed if x['rmsd_order_pair_00_delta']=='')['rmsd_order_pair_00_delta']='0.0'
            elif label=='changed_order_pair_01':next(x for x in changed if x['rmsd_usable_order_pair_bits']=='15')['rmsd_order_pair_01_delta']='999'
            elif label=='favorable_order_mean':next(x for x in changed if x['rmsd_usable_order_pair_bits']=='15')['rmsd_target_both_orders_mean']='1.5'
            elif label=='partial_envelope_promoted':next(x for x in changed if x['rmsd_usable_order_pair_bits']!='15')['rmsd_complete_order_envelope_min']='0.0'
            elif label=='changed_envelope_span':next(x for x in changed if x['rmsd_usable_order_pair_bits']=='15')['rmsd_complete_order_envelope_span']='0.0'
            elif label=='changed_pair_bitmap':changed[0]['rmsd_usable_order_pair_bits']='0'
            elif label=='changed_native_tm_dissimilarity':changed[0]['target_order0_native_tm_dissimilarity']='999'
            elif label=='changed_joint_quality':changed[0]['joint_mask_pass_bits']='63'
            elif label=='changed_same_model':changed[0]['target_same_model']='1'
            else:receipt['complete_envelopes']['rmsd:full']+=1
            zipped_table(data,changed);receipt['artifacts']['case_mask_contrasts.tsv.gz']=sha(data);write(rp,receipt);invoke(reader+[str(root/(label+'.json'))],False)
        data.write_bytes(saved);rp.write_bytes(original);rp.unlink();(out/'case_mask_contrasts.tsv.gz.tmp').write_text('interrupted')
        invoke(producer);invoke(reader+[str(root/'recovered.json')]);recovered=json.loads(rp.read_text())
        assert all(r[k]==recovered[k] for k in SUMMARY_FIELDS) and rows==list(csv.DictReader(gzip.open(data,'rt'),delimiter='\t'))
        return dict(status='passed_full_expanded_case_measurement_join_software_contracts',logical_cases=6,case_mask_rows=12,order_pair_cells=96,all_four_target_control_order_pairs_checked=True,both_physical_masks_and_quality_intersection_checked=True,distinct_gene_cases_on_shared_physical_pairs_retained=True,zero_and_null_outcomes_distinguished=True,incomplete_order_envelopes_not_promoted=True,independent_decimal_arithmetic_passed=True,completed_restart_refused=True,interrupted_full_replay_passed=True,rejected_rehashed_exports=labels,synthetic_prior_source_matching_measurement_and_journal_contracts=True,production_fits_or_journals_tested=False,scope=__doc__)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run();r['source_hashes']={n:sha(n) for n in ['scripts/full_expanded_case_measurement_sources.py','scripts/join_full_expanded_case_measurements.py','scripts/readback_full_expanded_case_measurements.py',__file__]}
    with a.output.open('x') as h:h.write(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
