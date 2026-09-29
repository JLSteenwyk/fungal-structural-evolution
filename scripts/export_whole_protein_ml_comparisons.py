#!/usr/bin/env python3
"""Export all audited whole-protein comparisons without inferential acceptance."""
import argparse
from collections import Counter
import fcntl
import gzip
import json
import math
from pathlib import Path
from whole_protein_support_registry import SupportRegistry
from screen_duplication_alignment_reuse import sha

PASS = 'ml_candidate_passed_numerical_optimization_checks'
REFINED = 'ml_refinement_passed_numerical_checks'


def comparison(left, right, relation):
    allowed = {'different_observations_no_direct_comparison', 'left_named_columns_nested_in_right',
               'right_named_columns_nested_in_left', 'same_observations_no_named_column_nesting', 'identical_named_design'}
    assert relation in allowed
    result = dict(log_likelihood_gain_right_vs_left=None, twice_log_likelihood_gain_right_vs_left=None,
                  nested_model_log_likelihood_gain=None, nested_added_fixed_coefficients=None,
                  numerical_likelihood_tolerance=None, nesting_status='not_assessed')
    if relation == 'different_observations_no_direct_comparison':
        result['comparison_status'] = relation
        return result
    assert left['records'] == right['records']
    if left['negative_profiled_ml'] is None or right['negative_profiled_ml'] is None:
        result['comparison_status'] = 'missing_fit_requires_review'
        return result
    a, b = left['negative_profiled_ml'], right['negative_profiled_ml']
    assert math.isfinite(a) and math.isfinite(b)
    gain = a-b; tolerance = 1e-7 + 1e-9*max(abs(a), abs(b))
    result.update(log_likelihood_gain_right_vs_left=gain, twice_log_likelihood_gain_right_vs_left=2*gain,
                  numerical_likelihood_tolerance=tolerance)
    violation = False
    if relation in ['left_named_columns_nested_in_right', 'right_named_columns_nested_in_left']:
        reduced, full = (left, right) if relation.startswith('left_') else (right, left)
        nested = reduced['negative_profiled_ml']-full['negative_profiled_ml']
        added = full['fixed_coefficients']-reduced['fixed_coefficients']; assert added > 0
        violation = nested < -tolerance
        result.update(nested_model_log_likelihood_gain=nested, nested_added_fixed_coefficients=added,
                      nesting_status='worse_full_model_requires_review' if violation else 'nondecreasing_within_numeric_tolerance')
    elif relation == 'identical_named_design':
        assert left['fixed_coefficients'] == right['fixed_coefficients']
        violation = abs(gain) > tolerance
        result['nesting_status'] = 'identical_design_likelihood_disagreement' if violation else 'identical_design_agreement_within_numeric_tolerance'
    else:
        result['nesting_status'] = 'no_named_column_nesting'
    result['comparison_status'] = ('nested_or_identical_likelihood_violation_requires_review' if violation else
        'optimization_review_required' if any(r['effective_status'] not in [PASS, REFINED] for r in [left, right]) else
        'zero_reference_support_review_required' if any(r['support_classification'] != 'zero_supported_to_numeric_tolerance' for r in [left, right]) else
        'numerically_checked_comparison_not_inferential_acceptance')
    return result


def selected_summary(entry, saved, registry, support):
    identifier, tree = entry['fit_input_id'], entry['tree']
    assert (saved['fit_input_id'], saved['tree']) == (identifier, tree)
    original = saved['payload']; spec = saved['specification']
    assert original['status'] == entry['status']
    if original['status'] != 'fit_error_requires_review':
        assert original['likelihood'] == 'ordinary_gaussian_ml'
    selected = original; kind = 'original_production'; selected_path = entry['path']; selected_sha = entry['sha256']
    effective = original['status']
    if (identifier, tree) in registry:
        recovery = registry[identifier, tree]
        assert recovery['original_production_fit'] == entry['path']
        assert recovery['original_production_sha256'] == entry['sha256']
        assert recovery['input_sha256'] == saved['input_sha256'] and recovery['original_status'] == original['status']
        candidate = recovery['candidate']; selected_path = candidate['path']; selected_sha = candidate['sha256']
        assert sha(selected_path) == selected_sha
        record = json.loads(Path(selected_path).read_text())
        assert (record['fit_input_id'], record['tree']) == (identifier, tree)
        kind = candidate['kind']
        if kind == 'audited_all_face_refinement':
            assert record['source_fit'] == entry['path'] and record['source_fit_sha256'] == entry['sha256']
            assert record['input_sha256'] == saved['input_sha256']
            selected = record['payload']; assert selected['status'] == REFINED
            assert selected['log1p_ratios'] == candidate['theta']
        else:
            assert kind == 'audited_full_face_start_recovery'
            assert record['status'] == 'full_face_start_recovery_passed_pending_readback' and all(record['checks'].values())
            assert record['selection'] == candidate['selection'] and record['selection']['theta'] == candidate['theta']
            selected = record['fitted']
        effective = REFINED
    p = len(spec['columns'])-1
    error = effective == 'fit_error_requires_review'
    if not error:
        assert len(selected['beta']) == p and selected['variance_profile_denominator'] == spec['records']
        assert math.isfinite(selected['negative_profiled_ml'])
    supported = support.resolve(identifier, saved['input_sha256'])
    return dict(fit_input_id=identifier, tree=tree, input_sha256=saved['input_sha256'],
                original_status=original['status'], effective_status=effective,
                original_negative_profiled_ml=original.get('negative_profiled_ml'),
                negative_profiled_ml=None if error else selected['negative_profiled_ml'],
                records=spec['records'], fixed_coefficients=p,
                original_production_fit=entry['path'], original_production_sha256=entry['sha256'],
                selected_source_kind=kind, selected_source_path=selected_path, selected_source_sha256=selected_sha,
                support_geometry_id=supported['geometry_id'], support_classification=supported['classification'],
                original_support_classification=supported['original_classification'], support_source_kind=supported['source_kind'])


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text()); ph = sha(args.plan)
    bindings = {str(args.plan):ph, **plan['pins']}
    root = Path(plan['fits']); fit_receipt = json.loads((root/'receipt.json').read_text())
    audit_root = Path(plan['fit_audit']); audit = json.loads((audit_root/'receipt.json').read_text())
    assert fit_receipt['status'] == 'complete_whole_protein_ml_dispositions_pending_full_audit'
    assert audit['status'] == 'complete_full_whole_protein_ml_output_audit'
    assert audit['source_receipt_sha256'] == sha(root/'receipt.json')
    assert fit_receipt['tree_fit_dispositions'] == audit['tree_fit_dispositions'] == 375350
    assert fit_receipt['plan_sha256'] == sha(plan['fit_plan'])
    for folder, receipt in [(root,fit_receipt),(audit_root,audit)]:
        bindings[str(folder/'receipt.json')] = sha(folder/'receipt.json')
        bindings.update({str(folder/name):h for name,h in receipt['artifacts'].items()})
    links = Path(plan['links']); link_receipt = json.loads((links/'receipt.json').read_text())
    link_proof = json.loads(Path(plan['link_proof']).read_text())
    assert link_proof['status'] == 'passed_full_whole_protein_comparison_input_link_readback'
    assert link_proof['source_receipt_sha256'] == sha(links/'receipt.json')
    bindings[str(links/'receipt.json')] = sha(links/'receipt.json')
    bindings.update({str(links/name):h for name,h in link_receipt['artifacts'].items()})
    recovery = json.loads(Path(plan['recoveries']).read_text())
    assert recovery['status'] == 'completed_frozen_five_fit_numerical_recovery_candidates'
    bindings.update(recovery['source_hashes'])
    registry = {(r['fit_input_id'],r['tree']):r for r in recovery['candidates']}
    assert len(registry) == recovery['fits'] == 5
    def verify():
        for path,digest in bindings.items(): assert sha(path) == digest,path
    verify()
    support = SupportRegistry(plan['support'], plan['support_proof'])
    bindings[str(Path(plan['support'])/'receipt.json')] = support.receipt_sha256
    bindings[plan['support_proof']] = support.proof_sha256
    output = Path(plan['output']); output.mkdir(parents=True,exist_ok=True)
    lock = (output/'run.lock').open('a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    stored = output/'run_plan.json'
    if stored.exists(): assert stored.read_bytes() == args.plan.read_bytes()
    else: stored.write_bytes(args.plan.read_bytes())
    audited = {}
    for line in (audit_root/'audit_manifest.jsonl').open():
        r = json.loads(line); key = r['fit_input_id'],r['tree']; assert key not in audited; audited[key] = r
    assert len(audited) == 375350
    summaries = {}; original_counts = Counter(); effective_counts = Counter(); summary_path = output/'fit_summaries.jsonl'
    temporary = output/'fit_summaries.jsonl.tmp'
    with temporary.open('w') as f:
        for line in (root/'fit_manifest.jsonl').open():
            entry = json.loads(line); key = entry['fit_input_id'],entry['tree']
            assert key not in summaries and sha(entry['path']) == entry['sha256']
            saved = json.loads(Path(entry['path']).read_text())
            assert saved['plan_sha256'] == fit_receipt['plan_sha256']
            assert audited[key]['status'] == entry['status']
            assert audited[key]['numerical_fit_verified'] == (entry['status'] != 'fit_error_requires_review')
            row = selected_summary(entry,saved,registry,support); summaries[key] = row
            original_counts[row['original_status']] += 1; effective_counts[row['effective_status']] += 1
            f.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
            bindings[entry['path']] = entry['sha256']
            if len(summaries)%10000 == 0: print('Whole-protein fit summaries',len(summaries),'/375350',flush=True)
    assert set(summaries) == set(audited) and dict(original_counts) == fit_receipt['status_counts']
    assert set(registry) <= set(summaries)
    temporary.replace(summary_path)
    trees = sorted({k[1] for k in summaries}); assert len(trees) == 5
    complete_keys = set(summaries)
    assert len({k[0] for k in complete_keys}) == 75070
    for identifier in {k[0] for k in complete_keys}: assert {(identifier,t) for t in trees} <= complete_keys
    counts = Counter(); relation_counts = Counter(); total = 0; artifacts = {'fit_summaries.jsonl':sha(summary_path),'run_plan.json':sha(stored)}
    summary_hash = sha(summary_path)
    for tree in trees:
        target = output/(tree+'.jsonl.gz'); checkpoint = output/(tree+'.receipt.json')
        if checkpoint.exists():
            record = json.loads(checkpoint.read_text())
            assert record['plan_sha256'] == ph and record['summary_sha256'] == summary_hash
            assert sha(target) == record['sha256'] and record['comparisons'] == 829440
        else:
            rows = 0; local = Counter(); relations = Counter(); temp = output/(tree+'.jsonl.gz.tmp')
            with gzip.open(temp,'wt',compresslevel=1) as f:
                for line in (links/'comparison_input_map.jsonl').open():
                    link = json.loads(line); left = summaries[link['left_input'],tree]; right = summaries[link['right_input'],tree]
                    row = dict(**link,tree=tree,left=left,right=right,**comparison(left,right,link['relation']))
                    f.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
                    rows += 1; local[row['comparison_status']] += 1; relations[link['relation']] += 1
            assert rows == 829440 and dict(relations) == link_receipt['relation_counts']
            temp.replace(target)
            record = dict(plan_sha256=ph,summary_sha256=summary_hash,sha256=sha(target),comparisons=rows,
                          status_counts=dict(local),relation_counts=dict(relations))
            checkpoint.write_text(json.dumps(record,indent=2)+'\n')
        counts.update(record['status_counts']); relation_counts.update(record['relation_counts']); total += record['comparisons']
        artifacts[target.name] = sha(target); artifacts[checkpoint.name] = sha(checkpoint)
        print('Exported full comparison tree',tree,total,'/4147200',flush=True)
    assert total == 4147200
    verify()
    result = dict(status='complete_full_whole_protein_ml_comparison_export_pending_readback',plan_sha256=ph,
                  scientific_eligibility=False,fit_summaries=len(summaries),comparisons=total,trees=trees,
                  original_fit_status_counts=dict(original_counts),effective_fit_status_counts=dict(effective_counts),
                  comparison_status_counts=dict(counts),relation_counts=dict(relation_counts),
                  audited_recovery_candidates_applied=len(registry),source_hashes=bindings,artifacts=artifacts,
                  scope='Full 829440 input comparisons across five trees. Original statuses and reviewed recovery choices retained. Different-observation comparisons have no likelihood gain; raw negative gains and review flags retained for comparable fits. No chi-square reference, p-values, AIC selection, global-optimum guarantee, calibrated uncertainty or scientific acceptance. Full serialized readback remains required.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__ == '__main__': main()
