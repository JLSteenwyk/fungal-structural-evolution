#!/usr/bin/env python3
"""Project all original-length screens onto full native contexts without reference reselection."""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from reference_context_coverage_sources import load_sources, POLICIES
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def side_coverage(work, mask, pairs, screens):
    status = work['measurement_disposition']; original = None
    if status == 'in_full_reference_measurement_design':
        original = pairs[work['pair_key'], mask]
        focal = ['a', 'b'][work['current_pair_focal_endpoint']]; ref = 'b' if focal == 'a' else 'a'
        assert (original['model_' + focal], int(original['version_' + focal])) == (work['duplicate_model'], work['duplicate_version'])
        assert (original['model_' + ref], int(original['version_' + ref])) == (work['reference_model'], work['reference_version'])
    result = dict(source_pair_mask_key=[work['pair_key'], mask] if original else None,
                  order_statuses=[original[f'order{o}_status'] for o in [0, 1]] if original else [None, None],
                  order_native_statuses=[original[f'order{o}_native_status'] for o in [0, 1]] if original else [None, None],
                  order_numerical_exclusions=[original[f'order{o}_numerical_exclusion_reasons'].split(';') if original[f'order{o}_numerical_exclusion_reasons'] else [] for o in [0, 1]] if original else [None, None],
                  screen_pass={}, screen_exclusions={})
    for spec in screens:
        name = spec['id']
        if original:
            assert original[name + '_pass'] in ['0', '1']
            passed = original[name + '_pass'] == '1'; why = original[name + '_exclusions'].split(';') if original[name + '_exclusions'] else []
            assert passed == (not why)
            if passed: assert all(s == 'aligned' for s in result['order_statuses'])
        else:
            assert status in ['reference_model_missing', 'identical_model_not_independent', 'outside_full_reference_measurement_design']
            passed = False; why = [status]
        result['screen_pass'][name] = passed; result['screen_exclusions'][name] = why
    return result


def project(source, pairs, screens):
    context = source['native_context']; eligible = context['parent_context_eligible']; guide = context['source_guide']
    result = dict(source_design=source, coverage_designs={})
    for design in ['availability', 'sequence_first']:
        refs = source['measurement_designs'][design]; coverage = []
        for ref in refs:
            sides = {side: {mask: side_coverage(ref['side_work'][side], mask, pairs, screens) for mask in ['full', 'plddt70']} for side in ['a', 'b']}
            coverage.append(dict(reference_gene=ref['reference']['reference_gene'], side_masks=sides))
        flags = {}
        for mask in ['full', 'plddt70']:
            flags[mask] = {}
            for spec in screens:
                name = spec['id']
                passing = [all(c['side_masks'][side][mask]['screen_pass'][name] for side in ['a', 'b']) for c in coverage]
                both = [passing[i] and all(v == 'both' for v in r['reference']['native_coorthology'].values()) for i, r in enumerate(refs)]
                lexical_pass = eligible and bool(refs) and passing[0]
                flags[mask][name] = [bool(lexical_pass), bool(lexical_pass and refs[0]['reference']['native_coorthology'][guide] == 'both'),
                                     bool(lexical_pass and both[0]), bool(eligible and any(both)), bool(eligible and bool(both) and all(both))]
        result['coverage_designs'][design] = dict(references=coverage, context_policy_flags=flags)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    contexts, table, upstream, coverage_receipt, bindings = load_sources(plan); bind(bindings, args.plan)
    pairs = {}
    with table.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['pair_key'], row['mask']; assert key not in pairs and key[1] in ['full', 'plddt70']; pairs[key] = row
    assert len(pairs) == 2 * plan['expected']['model_pairs']
    assert set(pairs) == {(key, mask) for key, _ in pairs for mask in ['full', 'plddt70']}
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    seen = set(); guides, baselines, policy_counts, side_counts = Counter(), Counter(), Counter(), Counter(); ties = logical = 0
    with contexts.open() as original, gzip.open(out / 'context_coverage.jsonl.gz', 'wt', compresslevel=1) as target:
        for line in original:
            source = json.loads(line); native = source['native_context']; key = native['source_guide'], native['source_row_number']
            assert key not in seen; seen.add(key); guides[key[0]] += 1
            result = project(source, pairs, plan['screens'])
            for design, d in result['coverage_designs'].items():
                refs = source['measurement_designs'][design]; parent = native['parent_context_eligible']; lexical = refs[0]['reference'] if refs else None
                base = key[0], design
                for name, flag in [('source_contexts', True), ('parent_eligible_contexts', parent), ('parent_eligible_lexical_gene_present', parent and bool(lexical)),
                                   ('parent_eligible_lexical_model_present', parent and bool(lexical and lexical['reference_model']))]: baselines[(*base, name)] += bool(flag)
                ties += len(refs); logical += 2 * len(refs)
                for ref in refs:
                    for work in ref['side_work'].values():
                        for mask in ['full', 'plddt70']: side_counts[f'{key[0]}|{design}|parent={int(parent)}|{mask}|{work["measurement_disposition"]}'] += 1
                for mask in ['full', 'plddt70']:
                    for spec in plan['screens']:
                        name = spec['id']; values = d['context_policy_flags'][mask][name]
                        assert len(values) == len(POLICIES) and all(isinstance(v, bool) for v in values)
                        assert not any(values) or parent
                        for policy, flag in zip(POLICIES, values): policy_counts[(*base, mask, name, policy)] += flag
            target.write(json.dumps(result, separators=(',', ':')) + '\n')
            if len(seen) % 25000 == 0: print('Full reference context coverage', len(seen), '/', plan['expected']['target_contexts'], flush=True)
    assert len(seen) == plan['expected']['target_contexts'] and dict(guides) == upstream['guide_contexts']
    assert seen == {(g, n) for g, count in guides.items() for n in range(1, count + 1)}
    assert ties == plan['expected']['reference_tie_records'] and logical == plan['expected']['duplicate_reference_links']
    assert sum(side_counts.values()) == 2 * logical
    for source_cell, count in upstream['measurement_disposition_counts'].items():
        guide, design, parent, status = source_cell.split('|')
        if status != 'no_reference_gene':
            for mask in ['full', 'plddt70']: assert side_counts[f'{guide}|{design}|{parent}|{mask}|{status}'] == count
    with (out / 'context_screen_counts.tsv').open('x') as handle:
        fields = ['guide', 'design', 'mask', 'screen', 'policy', 'source_contexts', 'parent_eligible_contexts', 'parent_eligible_lexical_gene_present', 'parent_eligible_lexical_model_present', 'passed_contexts']
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for key, count in sorted(policy_counts.items()):
            row = dict(zip(fields[:5], key), passed_contexts=count)
            for field in fields[5:9]: row[field] = baselines[key[0], key[1], field]
            writer.writerow(row)
    verify(bindings)
    summary = dict(target_contexts=len(seen), context_design_records=2 * len(seen), context_design_mask_records=4 * len(seen),
                   reference_tie_records=ties, duplicate_reference_links=logical, availability_side_links=upstream['availability_side_links_checked'],
                   context_screen_rows=24 * len(seen), context_policy_decisions=24 * len(seen) * len(POLICIES), side_screen_decisions=12 * logical,
                   summary_rows=len(policy_counts), guide_contexts=dict(guides), side_mask_work_counts=dict(side_counts))
    result = dict(status='complete_full_reference_context_coverage_pending_independent_readback', plan_sha256=sha(args.plan), **summary,
                  context_policy_flag_order=POLICIES, screens=plan['screens'], source_hashes=bindings,
                  artifacts={n: sha(out / n) for n in ['context_coverage.jsonl.gz', 'context_screen_counts.tsv']}, scientific_eligibility=False,
                  scope='Every source context/design/tie and all missing, identical and outside-design states retained with both masks and six coverage screens. Pair flags require the closed whole-reference order/coverage proof; numerical dispositions remain explicit and unmeasured order fields are null, not zero. Original parent eligibility gates all five lexical/any/all-tie policies, even for physical pairs measured elsewhere. Native own-guide/both-guide requirements separate; original lexical reference never changed. Every source field remains in source_design. Counts overlap across guides/designs/masks/screens/policies; these are pair-coverage/native-assignment eligibility diagnostics, not common-triad/domain/PAE/prediction/biological orthology/phylogeny or calibrated asymmetry acceptance.')
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts']}, indent=2), flush=True)


if __name__ == '__main__': main()
