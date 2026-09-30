#!/usr/bin/env python3
"""Independently reconstruct full sequence-first reference choice with leaf intervals and fsum paths."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from inventory_sequence_first_sister_references import sources
from readback_duplication_sister_references import index_tree, check_row, key
from screen_duplication_alignment_reuse import sha


def reconstruct_sequence_first(row, index, models):
    nodes, parents, spans, leaves = index
    name = row['gene_node']
    parent = parents[name]
    a, b = row['gene_a'], row['gene_b']
    if (set(leaves[slice(*spans[name])]) != {a, b} or len(nodes[name].clades) != 2
            or parents[a] != name or parents[b] != name):
        raise ValueError('Not a native terminal pair')
    paths = {}
    if parent is not None:
        lo, hi = spans[parent]
        da, db = spans[name]
        sister_leaves = leaves[lo:da] + leaves[db:hi]
        for gene in sister_leaves:
            if gene.split('_', 1)[0] == row['taxon_id']:
                continue
            edges = [nodes[name].branch_length]
            current = gene
            while current != parent:
                if current is None:
                    raise ValueError('Non-sister reference gene')
                edges.append(nodes[current].branch_length)
                current = parents[current]
            paths[gene] = math.fsum(edges)
    best = min(paths.values()) if paths else None
    genes = sorted(g for g in paths if math.isclose(paths[g], best, rel_tol=0, abs_tol=1e-12))
    chosen = genes[0] if genes else ''
    available_paths = [d for g, d in paths.items() if g in models]
    closest_modeled = min(available_paths) if available_paths else None
    assignments = []
    for gene in genes:
        model, version = models.get(gene, ('', ''))
        assignments.append(dict(gene=gene, model_id=model, version=version))
    if row['status'] not in {'provisional_reference_available', 'no_modeled_nonfocal_sister'}:
        status = 'ineligible_parent_context'
    elif not genes:
        status = 'no_nonfocal_sister_gene'
    elif chosen in models:
        status = 'sequence_first_lexical_reference_modeled'
    elif any(g in models for g in genes):
        status = 'lexical_reference_unmodeled_other_nearest_tie_modeled'
    elif row['chosen_reference_gene']:
        status = 'nearest_sequence_ties_unmodeled_farther_reference_used'
    else:
        status = 'no_modeled_sister_reference'
    return dict(sequence_first_nonfocal_genes=len(paths), sequence_first_nearest_genes=genes,
                sequence_first_tied_model_assignments=assignments, sequence_first_chosen_gene=chosen,
                sequence_first_chosen_model=models[chosen][0] if chosen in models else '',
                sequence_first_chosen_version=models[chosen][1] if chosen in models else '',
                sequence_first_minimum_distance='' if best is None else best,
                sequence_first_chosen_distance=paths[chosen] if chosen else '',
                nearest_modeled_minimum_distance='' if closest_modeled is None else closest_modeled,
                availability_distance_increment='' if closest_modeled is None else closest_modeled-best,
                sequence_first_ties=len(genes), sequence_first_modeled_ties=sum(g in models for g in genes),
                nearest_gene_set_overlap=len(set(genes) & set(json.loads(row['nearest_reference_genes']))),
                chosen_gene_matches_availability_choice=int(bool(chosen) and chosen == row['chosen_reference_gene']),
                sequence_first_status=status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins'],
                str(root / 'receipt.json'): sha(root / 'receipt.json')}
    for name, digest in receipt['artifacts'].items():
        bindings[str(root / name)] = digest

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source or artifact: ' + path)

    verify()
    source, original_receipt, reviewplan, models = sources(plan)
    if (receipt['status'] != 'complete_full_sequence_first_sister_reference_inventory_pending_readback'
            or receipt['plan_sha256'] != sha(args.plan)
            or receipt['original_reference_receipt_sha256'] != sha(source / 'receipt.json')
            or receipt['native_reference_readback_sha256'] != sha(plan['reference_readback'])
            or receipt['source_hashes'] != {str(args.plan): sha(args.plan), **plan['pins']}):
        raise ValueError('Unbound source receipt')
    reports = []
    for entry in reviewplan['guides']:
        guide = entry['guide']
        with (source / (guide + '_sister_references.tsv')).open() as handle:
            originals = list(csv.DictReader(handle, delimiter='\t'))
        with (root / (guide + '_sequence_first_references.tsv')).open() as handle:
            actual = list(csv.DictReader(handle, delimiter='\t'))
        observed = {key(row): row for row in actual}
        expected_keys = {key(row) for row in originals}
        if len(observed) != len(actual) or len(expected_keys) != len(originals) or set(observed) != expected_keys:
            raise ValueError('Wrong target universe')
        groups = defaultdict(list)
        for row in originals:
            record = observed[key(row)]
            if any(record[field] != value for field, value in row.items()):
                raise ValueError('Original source field differs')
            groups[row['family']].append(record)
        del originals, actual, observed
        seen, counts, contexts, agreement, matrix = set(), Counter(), Counter(), Counter(), Counter()
        delta = 0.0
        with Path(entry['trees']).open() as handle:
            for line in handle:
                family, newick = line.rstrip().split(': ', 1)
                if family not in groups:
                    continue
                if family in seen:
                    raise ValueError('Repeated native family tree')
                index = index_tree(newick)
                for row in groups[family]:
                    expected = reconstruct_sequence_first(row, index, models)
                    delta = max(delta, check_row(row, expected))
                    counts[expected['sequence_first_status']] += 1
                    contexts[row['status']] += 1
                    agreement[str(expected['chosen_gene_matches_availability_choice'])] += 1
                    matrix[row['status'] + '|' + expected['sequence_first_status']] += 1
                seen.add(family)
                if len(seen) % 1000 == 0:
                    print(guide, 'checked sequence-first families', len(seen), flush=True)
        if seen != set(groups):
            raise ValueError('Missing native family tree')
        expected_report = dict(guide=guide, candidates=sum(counts.values()), families=len(seen), counts=dict(counts),
                               source_context_counts=dict(contexts), chosen_gene_agreement=dict(agreement),
                               source_context_by_choice_status=dict(matrix))
        saved = next(r for r in receipt['guides'] if r['guide'] == guide)
        if saved != expected_report or saved['candidates'] != next(r['candidates'] for r in original_receipt['guides'] if r['guide'] == guide):
            raise ValueError('Wrong source disposition/summary counts')
        reports.append(dict(**expected_report, maximum_path_difference=delta))
        print(json.dumps(reports[-1]), flush=True)
    if {r['guide'] for r in reports} != {r['guide'] for r in receipt['guides']}:
        raise ValueError('Missing guide summary')
    verify()
    proof = dict(status='passed_full_sequence_first_sister_reference_native_readback',
                 producer_receipt_sha256=sha(root / 'receipt.json'), source_hashes=bindings, guides=reports,
                 scientific_eligibility=False,
                 scope='Complete original target universe and every original field preserved. Native leaf-interval subtraction and upward math.fsum paths independently reconstruct all nearest sequence-only ties, lexical choices, gene/model/version assignments, availability distance increments and complete context/choice matrices. Shares Newick parser/source-binding loader, not producer choice functions. No structural outcome, reference orthology, ancestor, calibrated uncertainty or asymmetry inference.')
    with args.output.open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    csv.field_size_limit(32*1024*1024)
    main()
