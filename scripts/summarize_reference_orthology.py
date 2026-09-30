#!/usr/bin/env python3
"""Report full premeasurement reference attrition from independently closed native assignments."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    completion = json.loads(Path(plan['completion']).read_text())
    assert completion['status'] == 'complete_verified_full_reference_native_orthology'
    source = Path(plan['source'])
    receipt_path = source / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert completion['source_hashes'][str(receipt_path)] == sha(receipt_path)
    assert completion['summary']['target_contexts'] == receipt['target_contexts'] == plan['target_contexts']
    assert len(completion['services']) == 2
    proof = json.loads(Path(plan['readback']).read_text())
    assert proof['status'] == 'passed_full_reference_native_orthology_readback'
    assert completion['source_hashes'][plan['readback']] == sha(plan['readback'])
    assert proof['producer_receipt_sha256'] == sha(receipt_path)
    context_path = source / 'contexts.jsonl'
    bindings.update({plan['completion']: sha(plan['completion']), plan['readback']: sha(plan['readback']),
                     str(receipt_path): sha(receipt_path), str(context_path): receipt['artifacts'][context_path.name]})
    assert completion['source_hashes'][str(context_path)] == bindings[str(context_path)]
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    totals, matrix = Counter(), Counter()
    contexts = 0
    for line in context_path.open():
        record = json.loads(line)
        guide, eligible = record['source_guide'], record['parent_context_eligible']
        contexts += 1
        for design in ['availability', 'sequence_first']:
            refs = record['designs'][design]
            lexical = refs[0] if refs else None
            if lexical:
                assert lexical['lexical_choice'] and all(not r['lexical_choice'] for r in refs[1:])
            model = bool(lexical and lexical['reference_model'])
            own = bool(lexical and lexical['native_coorthology'][guide] == 'both')
            both = bool(lexical and all(v == 'both' for v in lexical['native_coorthology'].values()))
            policies = dict(all_source_contexts=True, parent_context_eligible=eligible,
                            reference_gene_present=bool(refs), lexical_reference_model_present=model,
                            parent_eligible_reference_gene_present=eligible and bool(refs),
                            parent_eligible_lexical_model_present=eligible and model,
                            parent_eligible_lexical_native_own_coorthology=eligible and own,
                            parent_eligible_lexical_native_both_guide_coorthology=eligible and both,
                            parent_eligible_lexical_model_native_own_coorthology=eligible and model and own,
                            parent_eligible_lexical_model_native_both_guide_coorthology=eligible and model and both,
                            parent_eligible_all_ties_native_both_guide_coorthology=eligible and bool(refs) and
                                all(all(v == 'both' for v in r['native_coorthology'].values()) for r in refs),
                            parent_eligible_any_tie_native_both_guide_coorthology=eligible and
                                any(all(v == 'both' for v in r['native_coorthology'].values()) for r in refs))
            for policy, include in policies.items():
                totals[guide, design, policy] += int(include)
            matrix[guide, design, record['source']['status'], int(eligible),
                   int(bool(refs)), int(model),
                   lexical['native_coorthology']['profile'] if lexical else 'not_queried',
                   lexical['native_coorthology']['mafft'] if lexical else 'not_queried'] += 1
    assert contexts == plan['target_contexts']
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    with (out / 'context_policy_counts.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['guide', 'design', 'policy', 'contexts'])
        for key, value in sorted(totals.items()): writer.writerow([*key, value])
    with (out / 'lexical_reference_context_matrix.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['guide', 'design', 'source_context_status', 'parent_context_eligible',
                         'reference_gene_present', 'lexical_model_present',
                         'profile_native_coorthology', 'mafft_native_coorthology', 'contexts'])
        for key, value in sorted(matrix.items()): writer.writerow([*key, value])
    assert sum(matrix.values()) == 2 * contexts
    for guide, count in completion['summary']['guide_contexts'].items():
        for design in ['availability', 'sequence_first']:
            assert totals[guide, design, 'all_source_contexts'] == count
            assert sum(v for k, v in matrix.items() if k[:2] == (guide, design)) == count
            assert totals[guide, design, 'parent_eligible_all_ties_native_both_guide_coorthology'] <= totals[
                guide, design, 'parent_eligible_lexical_native_both_guide_coorthology'] <= totals[
                guide, design, 'parent_eligible_any_tie_native_both_guide_coorthology']
            assert totals[guide, design, 'parent_eligible_lexical_model_native_both_guide_coorthology'] <= totals[
                guide, design, 'parent_eligible_lexical_model_native_own_coorthology'] <= totals[
                guide, design, 'parent_eligible_lexical_model_present']
    verify()
    result = dict(status='complete_full_reference_native_assignment_attrition_summary',
                  plan_sha256=sha(args.plan), target_contexts=contexts, context_design_records=2 * contexts,
                  policy_rows=len(totals), matrix_rows=len(matrix),
                  source_hashes=bindings, artifacts={p.name: sha(p) for p in out.iterdir()}, scientific_eligibility=False,
                  scope='Complete premeasurement context counts and lexical model/native coorthology attrition, with '
                        'all parent, missing-gene/model and native disagreement states. Any/all tie policies are '
                        'reported separately; they do not replace the lexical choice or override parent exclusions. '
                        'Counts overlap across source guides, designs and nested policies, not independent events, '
                        'biological orthology validations, structural effects or calibrated tests.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
