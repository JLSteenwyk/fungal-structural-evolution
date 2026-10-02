#!/usr/bin/env python3
"""Qualify every full design under both uniform loading modes and five trees."""
from collections import Counter
import argparse
import csv
import fcntl
import gzip
import itertools
import json
from pathlib import Path
import shutil
import time
import numpy as np
from covariance_basis_context import ComponentKernelProducts
from full_covariance_qualification_sources import load,jsonl,cohort_rows,design_matrix,folded_operators,audit_id,covariance_status,SUMMARY_FIELDS,LINK_EXTRA
from full_entity_operator_sources import MODES
from full_expanded_model_input_sources import ORDERS
from full_expanded_model_design_sources import AXES,DEGREES,SETTING_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(path,stop_after_cohorts=None):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);root=Path(plan['output'])
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(),'Completed qualification cannot restart'
    state=dict(plan_sha256=sha(path),schema='expanded-uniform-covariance-qualification-v1');marker=root/'stage_plan.json'
    if marker.exists():assert json.loads(marker.read_text())==state
    else:marker.write_text(json.dumps(state,indent=2)+'\n')
    designs=jsonl(source['root']/'unique_designs.jsonl');statuses={};counts=Counter();audits=cohort_count=design_count=0
    with gzip.open(root/'design_covariance_audits.jsonl.gz','wt') as output:
        for cohort in source['cohorts']:
            rows=cohort_rows(source,cohort);contexts={};started=time.perf_counter()
            if len(rows):
                selected_factors={tree:source['factors'][tree][rows] for tree in plan['trees']}
                for mode in MODES:
                    bank=ComponentKernelProducts(source['labels'][rows],folded_operators(source,rows,mode),np.ones(len(rows)))
                    contexts[mode]={tree:bank.tree(selected_factors[tree]) for tree in plan['trees']}
            seen=set()
            for _ in range(30):
                design=next(designs);seen.add((design['order_contrast'],design['sequence_axis'],design['degree']))
                matrix=design_matrix(source,cohort,rows,design);design_count+=1
                for mode,tree in itertools.product(MODES,plan['trees']):
                    record=dict(audit_id=audit_id(source['contract'],design['design_id'],mode,tree),
                        source_contract=source['contract'],cohort_id=cohort['cohort_id'],design_id=design['design_id'],
                        loading_mode=mode,tree=tree,residual_diagonal='uniform_one',records=len(rows),
                        source_design_disposition=design['disposition'],folded_terms={
                            'target_node':'combined_with_uniform_residual',
                            'family':'zero_signed_endpoint_kernel' if mode=='signed' else 'combined_variance_family_intercept_plus_four_times_endpoint_family'},
                        active_column_indices=design['active_column_indices'],raw_design_sha256=design['raw_design_sha256'])
                    if design['disposition']!='full_rank_design':record.update(disposition=design['disposition'],numerical_audit=None)
                    else:
                        try:
                            value=contexts[mode][tree].audit(matrix)
                            record.update(disposition=covariance_status(value),numerical_audit=value)
                        except (ValueError,ArithmeticError) as error:
                            record.update(disposition='numerical_covariance_qualification_requires_review',numerical_audit=None,
                                error_type=type(error).__name__,error_message=str(error))
                    key=(design['design_id'],mode,tree);assert key not in statuses
                    statuses[key]=record['disposition'];counts[record['disposition']]+=1;audits+=1
                    output.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
            assert seen==set(itertools.product(ORDERS,AXES,DEGREES))
            cohort_count+=1
            print('full_covariance_qualified_cohort',cohort_count,'/',len(source['cohorts']),
                'records',len(rows),'seconds',time.perf_counter()-started,flush=True)
            if stop_after_cohorts==cohort_count:raise InterruptedError('Software interruption contract')
        assert next(designs,None) is None
    linkcounts=Counter();settings=links=0
    with gzip.open(source['root']/'model_settings.tsv.gz','rt') as f,gzip.open(root/'setting_audit_links.tsv.gz','wt') as g:
        reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==SETTING_FIELDS
        writer=csv.DictWriter(g,SETTING_FIELDS+LINK_EXTRA,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in reader:
            for mode,tree in itertools.product(MODES,plan['trees']):
                status=statuses[row['design_id'],mode,tree]
                combined=row['disposition'] if row['disposition']!='ready_for_working_covariance_fit' else status
                writer.writerow({**row,'loading_mode':mode,'tree':tree,'audit_id':audit_id(source['contract'],row['design_id'],mode,tree),
                    'covariance_disposition':status,'combined_disposition':combined})
                links+=1;linkcounts[combined]+=1
            settings+=1
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=settings,unique_cohorts=cohort_count,
        unique_designs=design_count,audit_rows=audits,setting_audit_links=links,
        audit_status_counts=dict(counts),link_status_counts=dict(linkcounts),trees=plan['trees'],loading_modes=MODES)
    assert settings==source['design_completion']['model_setting_rows']==plan['expected']['model_setting_rows']
    assert design_count==source['design_completion']['unique_designs'] and audits==design_count*len(MODES)*len(plan['trees'])
    assert links==settings*len(MODES)*len(plan['trees'])
    verify(bindings)
    artifacts={name:sha(root/name) for name in ['stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz']}
    receipt=dict(status='complete_full_uniform_covariance_qualification_pending_independent_readback',
        plan_sha256=sha(path),source_contract=source['contract'],**summary,artifacts=artifacts,source_hashes=bindings,
        scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
