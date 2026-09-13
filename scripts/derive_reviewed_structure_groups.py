#!/usr/bin/env python3
"""Retain only bidirectionally passing representative edges; defer failed members."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['clusters', 'review', 'annotations', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable output')
    receipts = {name: checked_receipt(getattr(args, name))
                for name in ['clusters', 'review', 'annotations']}
    if receipts['annotations']['cluster_receipt_sha256'] != sha(args.clusters / 'receipt.json'):
        raise ValueError('Annotation source mismatch')
    if receipts['review']['validation_receipt_sha256'] != receipts['annotations']['edge_validation_receipt_sha256']:
        raise ValueError('Review and annotations have different edge validation sources')
    members = read_table(args.clusters / 'model_cluster_membership.tsv')
    edges = read_table(args.review / 'edge_score_review.tsv')
    edge_by_pair = {(r['representative_input_id'], r['member_input_id']): r for r in edges}
    expected = {(r['representative_input_id'], r['cluster_input_id']) for r in members
                if r['representative_input_id'] != r['cluster_input_id']}
    if len(edge_by_pair) != len(edges) or set(edge_by_pair) != expected:
        raise ValueError('Incomplete or duplicate edge grid')
    identities = {r['cluster_input_id'] for r in members}
    if len(identities) != len(members) or len(members) != receipts['clusters']['models']:
        raise ValueError('Incomplete or duplicate model universe')
    retained, deferred, disposition = [], [], {}
    groups = defaultdict(list)
    for member in members:
        rep, model = member['representative_input_id'], member['cluster_input_id']
        if rep == model:
            status, retained_flag = 'representative_no_self_edge_test', True
        else:
            edge = edge_by_pair[rep, model]
            passes = edge['forward_exact_passes'] == 'True' and edge['reverse_exact_passes'] == 'True'
            if passes != (edge['exact_status'] == 'both_directions_pass'):
                raise ValueError('Inconsistent edge disposition')
            if passes and edge['any_exact_score_out_of_bounds'] != 'False':
                raise ValueError('Invalid score passed review')
            retained_flag = passes
            status = 'bidirectional_edge_pass' if passes else ('deferred_invalid_score' if edge['any_exact_score_out_of_bounds'] == 'True' else 'deferred_edge_criteria')
        row = dict(member, reviewed_group_id='RG_' + rep if retained_flag else '', disposition=status)
        disposition[model] = row
        (retained if retained_flag else deferred).append(row)
        if retained_flag:
            groups[rep].append(row)
    links = read_table(args.annotations / 'cluster_taxon_marker_links.tsv')
    if {r['cluster_input_id'] for r in links} != identities:
        raise ValueError('Annotation model universe mismatch')
    group_links = defaultdict(list)
    for row in links:
        assignment = disposition[row['cluster_input_id']]
        if row['representative_input_id'] != assignment['representative_input_id']:
            raise ValueError('Annotation representative mismatch')
        row.update(reviewed_group_id=assignment['reviewed_group_id'], disposition=assignment['disposition'])
        if row['reviewed_group_id']:
            group_links[row['representative_input_id']].append(row)
    summaries = []
    for rep, rows in sorted(groups.items()):
        if sum(r['cluster_input_id'] == rep for r in rows) != 1:
            raise ValueError('Missing representative')
        annotations = group_links[rep]
        summaries.append({'reviewed_group_id': 'RG_' + rep, 'representative_input_id': rep,
                          'models': len(rows), 'unique_sequences': len({r['sequence_sha256'] for r in rows}),
                          'taxa': len({r['taxon_id'] for r in annotations}),
                          'outgroup_taxa': len({r['taxon_id'] for r in annotations if r['study_role'] == 'outgroup'}),
                          'markers': ';'.join(sorted({r['marker'] for r in annotations})),
                          'sources': ';'.join(sorted({r['source'] for r in rows})),
                          'edge_scope': 'representative_star_only' if len(rows) > 1 else 'no_retained_nonself_edge',
                          'biological_family_status': 'not_established'})
    args.output.mkdir(parents=True)
    for name, rows in [('retained_membership.tsv', retained), ('deferred_members.tsv', deferred),
                       ('model_taxon_marker_links.tsv', links), ('group_summary.tsv', summaries)]:
        write_table(args.output / name, rows)
    result = {'status': 'complete_conservative_representative_group_derivation',
              'source_receipts': {name: sha(getattr(args, name) / 'receipt.json') for name in receipts},
              'script_sha256': sha(Path(__file__)), 'input_models': len(members),
              'retained_models': len(retained), 'deferred_models': len(deferred),
              'groups': len(groups), 'groups_without_nonself_edge': sum(len(v) == 1 for v in groups.values()),
              'dispositions': dict(Counter(r['disposition'] for r in disposition.values())),
              'model_taxon_marker_links': len(links),
              'interpretation': 'Conservative representative stars from fixed original alignments. Every nonself retained edge passes exact-option review in both directions. Deferred models have no derived group and are not novel singleton families. No member/member similarity guarantee, transitive equivalence, homology, orthology or function claim. Original memberships remain immutable. Source, confidence, domain and threshold sensitivities remain required.',
              'artifacts': {p.name: sha(p) for p in args.output.iterdir()}}
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
