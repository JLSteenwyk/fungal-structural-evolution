#!/usr/bin/env python3
"""Read every source summary/link and independently check comparison arithmetic."""
import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
from itertools import zip_longest
from export_whole_protein_ml_comparisons import selected_summary
from whole_protein_support_registry import SupportRegistry
from screen_duplication_alignment_reuse import sha


def verify_comparison(row, left, right, relation):
    assert row['left'] == left and row['right'] == right
    metrics = ['log_likelihood_gain_right_vs_left','twice_log_likelihood_gain_right_vs_left',
               'nested_model_log_likelihood_gain','nested_added_fixed_coefficients','numerical_likelihood_tolerance']
    if relation == 'different_observations_no_direct_comparison':
        assert all(row[k] is None for k in metrics)
        assert row['nesting_status'] == 'not_assessed' and row['comparison_status'] == relation
        return
    assert left['records'] == right['records']
    if left['negative_profiled_ml'] is None or right['negative_profiled_ml'] is None:
        assert all(row[k] is None for k in metrics)
        assert row['nesting_status'] == 'not_assessed' and row['comparison_status'] == 'missing_fit_requires_review'
        return
    a,b = left['negative_profiled_ml'],right['negative_profiled_ml']
    gain = math.fsum([a,-b]); tolerance = math.fsum([1e-7,1e-9*max(abs(a),abs(b))])
    assert row['log_likelihood_gain_right_vs_left'] == gain
    assert row['twice_log_likelihood_gain_right_vs_left'] == math.fsum([gain,gain])
    assert row['numerical_likelihood_tolerance'] == tolerance
    bad = False
    if relation == 'left_named_columns_nested_in_right':
        nested = gain; added = right['fixed_coefficients']-left['fixed_coefficients']
    elif relation == 'right_named_columns_nested_in_left':
        nested = -gain; added = left['fixed_coefficients']-right['fixed_coefficients']
    else:
        nested = added = None
    assert row['nested_model_log_likelihood_gain'] == nested and row['nested_added_fixed_coefficients'] == added
    if nested is not None:
        assert added > 0
        bad = nested < -tolerance
        nesting = 'worse_full_model_requires_review' if bad else 'nondecreasing_within_numeric_tolerance'
    elif relation == 'identical_named_design':
        assert left['fixed_coefficients'] == right['fixed_coefficients']
        bad = abs(gain) > tolerance
        nesting = 'identical_design_likelihood_disagreement' if bad else 'identical_design_agreement_within_numeric_tolerance'
    else:
        assert relation == 'same_observations_no_named_column_nesting'
        nesting = 'no_named_column_nesting'
    assert row['nesting_status'] == nesting
    if bad: status = 'nested_or_identical_likelihood_violation_requires_review'
    elif {left['effective_status'],right['effective_status']} - {'ml_candidate_passed_numerical_optimization_checks','ml_refinement_passed_numerical_checks'}:
        status = 'optimization_review_required'
    elif {left['support_classification'],right['support_classification']} != {'zero_supported_to_numeric_tolerance'}:
        status = 'zero_reference_support_review_required'
    else: status = 'numerically_checked_comparison_not_inferential_acceptance'
    assert row['comparison_status'] == status


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',required=True,type=Path)
    args=parser.parse_args();config=json.loads(args.plan.read_text())
    plan=json.loads(Path(config['source_plan']).read_text());root=Path(plan['output'])
    rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_whole_protein_ml_comparison_export_pending_readback'
    assert receipt['scientific_eligibility'] is False and receipt['plan_sha256']==sha(config['source_plan'])
    bindings={str(args.plan):sha(args.plan),**config['pins'],**receipt['source_hashes'],str(rp):sha(rp)}
    bindings.update({str(root/name):h for name,h in receipt['artifacts'].items()})
    def verify():
        for path,digest in bindings.items():assert sha(path)==digest,path
    verify()
    support=SupportRegistry(plan['support'],plan['support_proof'])
    recovery=json.loads(Path(plan['recoveries']).read_text());registry={(r['fit_input_id'],r['tree']):r for r in recovery['candidates']}
    originals={}
    for line in (Path(plan['fits'])/'fit_manifest.jsonl').open():
        r=json.loads(line);key=r['fit_input_id'],r['tree'];assert key not in originals;originals[key]=r
    summaries={};original_counts=Counter();effective_counts=Counter()
    for line in (root/'fit_summaries.jsonl').open():
        saved=json.loads(line);key=saved['fit_input_id'],saved['tree'];assert key not in summaries and key in originals
        source=originals[key];assert sha(source['path'])==source['sha256']
        source_payload=json.loads(Path(source['path']).read_text())
        expected=selected_summary(source,source_payload,registry,support)
        assert saved==expected
        summaries[key]=saved;original_counts[saved['original_status']]+=1;effective_counts[saved['effective_status']]+=1
    assert set(summaries)==set(originals) and len(summaries)==receipt['fit_summaries']==375350
    assert dict(original_counts)==receipt['original_fit_status_counts'] and dict(effective_counts)==receipt['effective_fit_status_counts']
    total=0;counts=Counter();relations=Counter();used=set()
    trees=sorted({key[1] for key in summaries});assert trees==receipt['trees'] and len(trees)==5
    for tree in trees:
        local=Counter();relation_counts=Counter();number=0
        with gzip.open(root/(tree+'.jsonl.gz'),'rt') as exported,(Path(plan['links'])/'comparison_input_map.jsonl').open() as links:
            for source_line,output_line in zip_longest(links,exported):
                assert source_line is not None and output_line is not None
                link=json.loads(source_line);row=json.loads(output_line)
                assert {k:row[k] for k in link}==link and row['tree']==tree
                assert set(row)==set(link)|{'tree','left','right','log_likelihood_gain_right_vs_left','twice_log_likelihood_gain_right_vs_left','nested_model_log_likelihood_gain','nested_added_fixed_coefficients','numerical_likelihood_tolerance','nesting_status','comparison_status'}
                keys=[(link[name],tree) for name in ['left_input','right_input']]
                left,right=[summaries[k] for k in keys];used.update(keys)
                verify_comparison(row,left,right,link['relation'])
                local[row['comparison_status']]+=1;relation_counts[link['relation']]+=1;number+=1
        tree_receipt=json.loads((root/(tree+'.receipt.json')).read_text())
        assert number==tree_receipt['comparisons']==829440
        assert dict(local)==tree_receipt['status_counts'] and dict(relation_counts)==tree_receipt['relation_counts']
        counts.update(local);relations.update(relation_counts);total+=number
        print('Read back whole-protein comparison tree',tree,total,'/4147200',flush=True)
    assert total==receipt['comparisons']==4147200 and used==set(summaries)
    assert dict(counts)==receipt['comparison_status_counts'] and dict(relations)==receipt['relation_counts']
    verify()
    result=dict(status='passed_full_whole_protein_ml_comparison_export_readback',comparisons=total,fit_summaries=len(summaries),
                trees=trees,comparison_status_counts=dict(counts),relation_counts=dict(relations),
                source_receipt_sha256=sha(rp),source_plan_sha256=sha(config['source_plan']),checker_sha256=sha(__file__),
                scope='Every source fit summary, recovery choice and setting/tree comparison read back. Separate arithmetic/orientation/nesting/decision implementation; selected-summary reader and support gate are shared. Different observations, errors, unresolved flags and negative nested gains retained. No calibrated likelihood-ratio distribution, p-values, biological significance or global-optimum claim.')
    with Path(config['output']).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
