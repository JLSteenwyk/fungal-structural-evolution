#!/usr/bin/env python3
"""Reconstruct every contracted/pruned input from audited source split sets."""
import argparse
from collections import Counter
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path

import dendropy
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha
from species_coalescent_sources import load


def parse(text, index, clean=False):
    tree = dendropy.Tree.get(data=text, schema='newick', rooting='force-unrooted', preserve_underscores=True)
    tips = [tip.taxon.label for tip in tree.leaf_node_iter()]
    assert len(tips) == len(set(tips)) >= 4 and set(tips) <= set(index)
    universe = sum(1 << index[tip] for tip in tips)
    desc, splits = {}, {}
    for node in tree.postorder_node_iter():
        if clean:
            assert node.edge_length is None and node.label is None and not node.annotations and not node.edge.annotations
        desc[node] = 1 << index[node.taxon.label] if node.is_leaf() else sum(desc[child] for child in node.child_node_iter())
        if not node.is_leaf(): assert len(node.child_nodes()) >= 2
        if node.parent_node is None: continue
        mask = min(desc[node], universe ^ desc[node])
        if mask.bit_count() > 1 and (universe ^ mask).bit_count() > 1:
            label = Decimal(node.label) if node.label is not None else None
            if mask in splits:
                # An arbitrary degree-two Newick root can represent both sides
                # of the same unrooted edge; clean inputs have no support labels.
                assert clean and len(tree.seed_node.child_nodes()) == 2 and node.parent_node is tree.seed_node
            splits[mask] = label
    return tips, universe, splits


def expected(source_rows, retained_mask, cutoff, index):
    counts = Counter()
    surviving = set()
    for row in source_rows:
        counts['original_internal'] += 1
        value = row['support']
        if cutoff is not None and (value is None or value < cutoff):
            counts['contracted_missing_support' if value is None else 'contracted_below_support'] += 1
            continue
        counts['eligible_original_internal'] += 1
        mask = sum(1 << index[tip] for tip in row['side']) & retained_mask
        key = min(mask, retained_mask ^ mask)
        if key.bit_count() > 1 and (retained_mask ^ key).bit_count() > 1: surviving.add(key)
    return surviving, {key: counts[key] for key in ['original_internal', 'contracted_missing_support',
                                                  'contracted_below_support', 'eligible_original_internal']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    universe, markers, cohorts, support, bindings = load(plan, args.plan)
    index = {taxon: i for i, taxon in enumerate(sorted(universe))}
    originals = {}
    for alignment, items in markers.items():
        for item in items:
            tips, mask, splits = parse(Path(item['path']).read_text(), index)
            assert len(tips) == item['taxa'] and len(splits) == item['taxa'] - 3
            audited = {}
            for row in support[alignment][item['marker']]:
                side = sum(1 << index[tip] for tip in row['side'])
                key = min(side, mask ^ side)
                assert key not in audited
                audited[key] = row['support']
            assert audited == splits
            originals[alignment, item['marker']] = dict(tips=set(tips), mask=mask)
    root = Path(plan['output'])
    rp = root / 'receipt.json'
    producer = json.loads(rp.read_text())
    assert producer['status'] == 'complete_full_coalescent_input_preparation_pending_independent_readback'
    assert producer['plan_sha256'] == sha(args.plan) and producer['scientific_eligibility'] is False
    assert producer['cases'] == 30 and producer['marker_states'] == 3750 and producer['original_gene_split_support_decisions'] == 1781685
    for p, d in producer['source_hashes'].items(): bind(bindings, p, d)
    bind(bindings, rp)
    for name, digest in producer['artifacts'].items(): bind(bindings, root / name, digest)
    config = json.loads((root / 'config.json').read_text())
    assert config['plan_sha256'] == sha(args.plan)
    for p, d in config['source_hashes'].items(): bind(bindings, p, d)
    ledger = {}
    with gzip.open(root / 'marker_cohort_support_inputs.jsonl.gz', 'rt') as handle:
        for line in handle:
            row = json.loads(line)
            key = row['case'], row['marker']
            assert key not in ledger
            ledger[key] = row
    assert len(ledger) == 3750
    expected_grid = {(alignment, cohort, policy) for alignment in markers for cohort in cohorts for policy in plan['support_policies']}
    grid = {(row['alignment'], row['cohort'], row['support_policy']): row for row in producer['summaries']}
    assert len(grid) == len(producer['summaries']) == 30 and set(grid) == expected_grid
    states = decisions = split_states = 0
    case_summaries = []
    for (alignment, cohort, policy), summary in grid.items():
        name = alignment + '-' + cohort + '-' + policy
        spec = cohorts[cohort]
        folder = root / name
        roles = dict(taxa=sorted(spec['taxa']), outgroups=sorted(spec['outgroups']), roles=spec['roles'])
        assert summary == dict(case=name, alignment=alignment, cohort=cohort, support_policy=policy,
                               cutoff=plan['support_policies'][policy], markers=125, expected_taxa=len(spec['taxa']), roles=spec['roles'],
                               genes=str(folder / 'genes.tree'), gene_order=str(folder / 'gene_order.json'), taxa_roles=str(folder / 'taxa_roles.json'))
        assert json.loads((folder / 'taxa_roles.json').read_text()) == roles
        assert json.loads((folder / 'gene_order.json').read_text()) == [item['marker'] for item in markers[alignment]]
        lines = (folder / 'genes.tree').read_text().splitlines()
        assert len(lines) == 125 and all(line.endswith(';') and line.count(';') == 1 for line in lines)
        totals = Counter()
        seen = set()
        for item, line in zip(markers[alignment], lines):
            origin = originals[alignment, item['marker']]
            intended = origin['tips'] & spec['taxa']
            tips, mask, splits = parse(line, index, clean=True)
            assert set(tips) == intended
            wanted, counts = expected(support[alignment][item['marker']], mask, plan['support_policies'][policy], index)
            assert set(splits) == wanted, (name, item['marker'])
            row = ledger.pop((name, item['marker']))
            assert row == dict(case=name, alignment=alignment, cohort=cohort, support_policy=policy,
                               cutoff=plan['support_policies'][policy], marker=item['marker'], source=item['path'], source_tree_sha256=item['tree_sha256'],
                               source_taxa=item['taxa'], retained_taxa=len(tips), retained_tips=sorted(tips),
                               gene_tree_sha256=hashlib.sha256(line.encode()).hexdigest(), **counts)
            states += 1
            decisions += counts['original_internal']
            split_states += len(splits)
            totals.update(counts)
            seen.update(tips)
        assert seen == spec['taxa']
        case_summaries.append(dict(case=name, markers=125, taxa=len(seen), roles=spec['roles'], original_gene_split_support_decisions=totals['original_internal'],
                                   contracted_missing_support=totals['contracted_missing_support'], contracted_below_support=totals['contracted_below_support']))
        print('independent_coalescent_input_case_verified', name, len(case_summaries), '/30', flush=True)
    assert not ledger and states == 3750 and decisions == 1781685
    verify(bindings)
    result = dict(status='passed_full_coalescent_input_independent_split_projection_readback',
                  plan_sha256=sha(args.plan), producer_receipt_sha256=sha(rp), cases=30, alignments=2, cohorts=5, support_policies=3,
                  marker_states=states, original_gene_split_support_decisions=decisions, retained_internal_split_states=split_states,
                  case_summaries=case_summaries, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'cases', 'marker_states', 'original_gene_split_support_decisions', 'retained_internal_split_states']}), flush=True)


if __name__ == '__main__':
    main()
