#!/usr/bin/env python3
"""Compare every ML/consensus split across all four crossed alignment/guide runs."""
import argparse
import csv
import itertools
import json
from pathlib import Path
from Bio import Phylo
from four_run_pmsf_sources import load_sources, PAIR_FIELDS, PRESENCE_FIELDS, CONFLICT_FIELDS, BOUNDARY_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def canonical(side, universe):
    return min(tuple(sorted(side)), tuple(sorted(universe - set(side))), key=lambda x: (len(x), x))


def rule(kind):
    return 'SH_aLRT80_and_empirical_UFB95' if kind == 'ml' else 'empirical_UFB95_SH_aLRT_unavailable'


def supported(row, kind):
    return float(row['empirical_ufboot_percent']) >= 95 and (kind == 'consensus' or float(row['sh_alrt_percent']) >= 80)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', required=True, type=Path)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); sources, manifest, bindings = load_sources(plan, args.plan)
    universe = {r['taxon_id'] for r in manifest}; views = {}; kinds = {}; runs = {}
    for label, source in sources.items():
        for kind, file in [('ml', 'pmsf.treefile'), ('consensus', 'pmsf.contree')]:
            name = label + ':' + kind; tree = Phylo.read(Path(source['spec']['run']) / file, 'newick')
            assert len(tree.get_terminals()) == len(universe) and {n.name for n in tree.get_terminals()} == universe
            split_rows = {}; raw = {}
            for node in tree.get_nonterminals():
                side = {n.name for n in node.get_terminals()}
                if len(side) in [0, 1, len(universe) - 1, len(universe)]: continue
                key = canonical(side, universe); assert key not in raw
                raw[key] = node
            for row in source['rows']:
                if row['tree'] != kind: continue
                key = canonical(set(json.loads(row['split_taxa_json'])), universe); assert key not in split_rows and key in raw
                assert float(row['branch_length']) == raw[key].branch_length
                split_rows[key] = row
            assert len(split_rows) == len(raw) == plan['expected']['internal_splits_per_view']
            views[name] = split_rows; kinds[name] = kind; runs[name] = label
    union = sorted(set().union(*(set(v) for v in views.values())))
    presence = []
    for key in union:
        for name, splits in views.items():
            row = splits.get(key)
            presence.append(dict(view=name, run=runs[name], tree_type=kinds[name], split_taxa_json=json.dumps(key), present=row is not None,
                                 branch_length=row['branch_length'] if row else '', sh_alrt=row['sh_alrt_percent'] if row else '', empirical_ufb=row['empirical_ufboot_percent'] if row else ''))
    pair_grid = [(a, b) for a, b in itertools.combinations(views, 2) if kinds[a] == kinds[b] or runs[a] == runs[b]]
    assert len(pair_grid) == 16
    comparison = []; conflicts = []
    for a, b in pair_grid:
        left, right = views[a], views[b]; common = set(left) & set(right); nconf = high_count = 0
        scope = 'within_run_ML_vs_consensus' if runs[a] == runs[b] else ('cross_run_' + kinds[a])
        for x, y in itertools.product(sorted(set(left) - common), sorted(set(right) - common)):
            xs, ys = set(x), set(y); cells = [xs & ys, xs - ys, ys - xs, universe - (xs | ys)]
            if not all(cells): continue
            nconf += 1; sx, sy = left[x], right[y]; high = supported(sx, kinds[a]) and supported(sy, kinds[b]); high_count += high
            conflicts.append(dict(view_a=a, view_b=b, split_a_taxa_json=json.dumps(x), split_b_taxa_json=json.dumps(y),
                                  support_rule_a=rule(kinds[a]), support_rule_b=rule(kinds[b]), sh_alrt_a=sx['sh_alrt_percent'], sh_alrt_b=sy['sh_alrt_percent'],
                                  empirical_ufb_a=sx['empirical_ufboot_percent'], empirical_ufb_b=sy['empirical_ufboot_percent'],
                                  both_support_criteria_met=high, witness_quartet=';'.join(min(cell) for cell in cells)))
        distance = len(left) + len(right) - 2 * len(common)
        comparison.append(dict(view_a=a, view_b=b, comparison_scope=scope, support_rule_a=rule(kinds[a]), support_rule_b=rule(kinds[b]),
                               shared_internal_splits=len(common), unique_internal_splits_a=len(left) - len(common), unique_internal_splits_b=len(right) - len(common),
                               rf_distance=distance, normalized_rf=distance / (len(left) + len(right)), incompatible_split_pairs=nconf, both_supported_incompatible_pairs=high_count))
    outgroups = {r['taxon_id'] for r in manifest if r['study_role'] == 'outgroup'}; boundary_key = canonical(outgroups, universe)
    boundary = []
    for name, splits in views.items():
        row = splits.get(boundary_key)
        boundary.append(dict(view=name, ingroup_count=501, outgroup_count=25, boundary_split_present=row is not None,
                             sh_alrt=row['sh_alrt_percent'] if row else '', empirical_ufb=row['empirical_ufboot_percent'] if row else '', boundary_taxa_json=json.dumps(boundary_key)))
    out = Path(plan['output']); out.mkdir(exist_ok=False, parents=True)
    for name, records, fields in [('comparisons.tsv', comparison, PAIR_FIELDS), ('split_presence.tsv', presence, PRESENCE_FIELDS), ('conflicts.tsv', conflicts, CONFLICT_FIELDS), ('rooting_boundary.tsv', boundary, BOUNDARY_FIELDS)]:
        with (out / name).open('x') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(records)
    verify(bindings)
    result = dict(status='complete_full_four_run_pmsf_ML_and_consensus_sensitivity_pending_readback', plan_sha256=sha(args.plan), taxa=len(universe), runs=4, tree_views=8,
                  internal_splits_per_view=523, comparison_rows=len(comparison), split_presence_rows=len(presence), incompatible_pairs=len(conflicts),
                  shared_all_ml=len(set.intersection(*(set(v) for name, v in views.items() if kinds[name] == 'ml'))),
                  shared_all_consensus=len(set.intersection(*(set(v) for name, v in views.items() if kinds[name] == 'consensus'))),
                  shared_all_eight=len(set.intersection(*(set(v) for v in views.values()))), boundary_views_with_role_split=sum(r['boundary_split_present'] for r in boundary),
                  comparisons=comparison, source_hashes=bindings, artifacts={name: sha(out / name) for name in ['comparisons.tsv', 'split_presence.tsv', 'conflicts.tsv', 'rooting_boundary.tsv']},
                  scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'artifacts', 'scope']}, indent=2), flush=True)


if __name__ == '__main__': main()
