#!/usr/bin/env python3
"""Prepare every marker under five cohorts and three gene-support policies.

SH-aLRT thresholds are sensitivity settings, not bootstrap thresholds. Missing
support is separately recorded and contracted only in support-filtered inputs.
Contraction precedes taxon pruning. Gene input branch lengths and labels are
omitted: ASTRAL uses unrooted topologies, not sequence branch lengths.
"""
import argparse
from collections import Counter
import copy
import gzip
import json
from pathlib import Path
import re

from Bio import Phylo
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha
from species_coalescent_sources import load


def topology(tree):
    def write(node):
        if node.is_terminal():
            assert re.fullmatch(r'[A-Za-z0-9_]+', node.name)
            return node.name
        return '(' + ','.join(write(child) for child in node.clades) + ')'
    return write(tree.root) + ';'


def transform(tree, retained, cutoff):
    result = copy.deepcopy(tree)
    counts = Counter()
    for node in list(result.get_nonterminals(order='postorder')):
        if node is result.root:
            continue
        counts['original_internal'] += 1
        if cutoff is not None and (node.confidence is None or node.confidence < cutoff):
            reason = 'contracted_missing_support' if node.confidence is None else 'contracted_below_support'
            counts[reason] += 1
            result.collapse(node)
        else:
            counts['eligible_original_internal'] += 1
    for name in sorted({tip.name for tip in result.get_terminals()} - retained):
        result.prune(next(tip for tip in result.get_terminals() if tip.name == name))
    names = [tip.name for tip in result.get_terminals()]
    assert len(names) == len(set(names)) >= 4 and set(names) <= retained
    return result, {key: counts[key] for key in ['original_internal', 'contracted_missing_support',
                                               'contracted_below_support', 'eligible_original_internal']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    universe, markers, cohorts, support, bindings = load(plan, args.plan)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    (out / 'config.json').write_text(json.dumps(dict(plan_sha256=sha(args.plan), source_hashes=bindings), indent=2) + '\n')
    marker_lookup = {}
    for alignment, items in markers.items():
        marker_lookup[alignment] = {}
        for item in items:
            tree = Phylo.read(item['path'], 'newick')
            assert len(tree.get_terminals()) == item['taxa']
            marker_lookup[alignment][item['marker']] = tree
    summaries = []
    total_states = support_decisions = 0
    with gzip.open(out / 'marker_cohort_support_inputs.jsonl.gz', 'xt') as ledger:
        for alignment, items in markers.items():
            for cohort, spec in cohorts.items():
                for policy, cutoff in plan['support_policies'].items():
                    name = alignment + '-' + cohort + '-' + policy
                    folder = out / name
                    folder.mkdir()
                    order = []
                    observed_union = set()
                    with (folder / 'genes.tree').open('x') as handle:
                        for item in items:
                            tree, counts = transform(marker_lookup[alignment][item['marker']], spec['taxa'], cutoff)
                            text = topology(tree)
                            handle.write(text + '\n')
                            tips = sorted(tip.name for tip in tree.get_terminals())
                            observed_union.update(tips)
                            row = dict(case=name, alignment=alignment, cohort=cohort, support_policy=policy,
                                       cutoff=cutoff, marker=item['marker'], source=item['path'], source_tree_sha256=item['tree_sha256'],
                                       source_taxa=item['taxa'], retained_taxa=len(tips), retained_tips=tips,
                                       gene_tree_sha256=__import__('hashlib').sha256(text.encode()).hexdigest(), **counts)
                            ledger.write(json.dumps(row, separators=(',', ':')) + '\n')
                            order.append(item['marker'])
                            total_states += 1
                            support_decisions += counts['original_internal']
                    assert observed_union == spec['taxa'], (name, sorted(spec['taxa'] - observed_union))
                    (folder / 'gene_order.json').write_text(json.dumps(order, indent=2) + '\n')
                    (folder / 'taxa_roles.json').write_text(json.dumps(dict(taxa=sorted(spec['taxa']), outgroups=sorted(spec['outgroups']), roles=spec['roles']), indent=2) + '\n')
                    record = dict(case=name, alignment=alignment, cohort=cohort, support_policy=policy,
                                  cutoff=cutoff, markers=125, expected_taxa=len(spec['taxa']), roles=spec['roles'],
                                  genes=str(folder / 'genes.tree'), gene_order=str(folder / 'gene_order.json'),
                                  taxa_roles=str(folder / 'taxa_roles.json'))
                    summaries.append(record)
                    print('full_coalescent_input_case_prepared', name, len(summaries), '/30', flush=True)
    assert len(summaries) == 30 and total_states == 3750 and support_decisions == 1781685
    artifacts = {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}
    verify(bindings)
    result = dict(status='complete_full_coalescent_input_preparation_pending_independent_readback',
                  plan_sha256=sha(args.plan), cases=30, alignments=2, cohorts=5, support_policies=3,
                  marker_states=3750, original_gene_split_support_decisions=support_decisions,
                  summaries=summaries, artifacts=artifacts, source_hashes=bindings,
                  scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'cases', 'marker_states', 'original_gene_split_support_decisions']}), flush=True)


if __name__ == '__main__':
    main()
