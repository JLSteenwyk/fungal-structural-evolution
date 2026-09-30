#!/usr/bin/env python3
"""Choose nearest sister genes before checking model availability for the full target cohort."""
import argparse
import csv
import json
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from inventory_duplication_sister_references import tree_index
from screen_duplication_alignment_reuse import sha


def classify(row, ties, assignments):
    # Parent eligibility is fixed by the independently checked native source.
    if row['status'] not in ['provisional_reference_available', 'no_modeled_nonfocal_sister']:
        return 'ineligible_parent_context'
    if not ties:
        return 'no_nonfocal_sister_gene'
    if assignments[0]['model_id']:
        return 'sequence_first_lexical_reference_modeled'
    if any(r['model_id'] for r in assignments):
        return 'lexical_reference_unmodeled_other_nearest_tie_modeled'
    if row['chosen_reference_gene']:
        return 'nearest_sequence_ties_unmodeled_farther_reference_used'
    return 'no_modeled_sister_reference'


def sequence_first(row, index, models):
    nodes, parents, _ = index
    node = nodes[row['gene_node']]
    if len(node.clades) != 2 or sorted(c.name for c in node.clades if not c.clades) != sorted([row['gene_a'], row['gene_b']]):
        raise ValueError('Not an exact native terminal pair')
    parent = parents[node.name]
    distances = {}
    if parent is not None:
        stack = [(c, c.branch_length) for c in nodes[parent].clades if c is not node]
        while stack:
            tip, path = stack.pop()
            if tip.clades:
                stack.extend((child, path + child.branch_length) for child in tip.clades)
            elif tip.name.split('_', 1)[0] != row['taxon_id']:
                distances[tip.name] = node.branch_length + path
    minimum = min(distances.values(), default=None)
    ties = sorted(g for g, d in distances.items() if abs(d-minimum) <= 1e-12)
    assignments = [dict(gene=g, model_id=models[g][0] if g in models else '',
                        version=models[g][1] if g in models else '') for g in ties]
    chosen = ties[0] if ties else ''
    available = [distances[g] for g in distances if g in models]
    minimum_modeled = min(available, default=None)
    observed_ties = set(json.loads(row['nearest_reference_genes']))
    common = len(set(ties) & observed_ties)
    return dict(sequence_first_nonfocal_genes=len(distances), sequence_first_nearest_genes=ties,
                sequence_first_tied_model_assignments=assignments, sequence_first_chosen_gene=chosen,
                sequence_first_chosen_model=assignments[0]['model_id'] if assignments else '',
                sequence_first_chosen_version=assignments[0]['version'] if assignments else '',
                sequence_first_minimum_distance='' if minimum is None else minimum,
                sequence_first_chosen_distance='' if not chosen else distances[chosen],
                nearest_modeled_minimum_distance='' if minimum_modeled is None else minimum_modeled,
                availability_distance_increment='' if minimum_modeled is None else minimum_modeled-minimum,
                sequence_first_ties=len(ties), sequence_first_modeled_ties=sum(bool(r['model_id']) for r in assignments),
                nearest_gene_set_overlap=common,
                chosen_gene_matches_availability_choice=int(bool(chosen) and chosen == row['chosen_reference_gene']),
                sequence_first_status=classify(row, ties, assignments))


def sources(plan):
    sourceplan = json.loads(Path(plan['reference_plan']).read_text())
    reference = Path(plan['references'])
    receipt = json.loads((reference / 'receipt.json').read_text())
    proof = json.loads(Path(plan['reference_readback']).read_text())
    if (receipt['status'] != 'complete_duplication_sister_reference_inventory'
            or receipt['plan_sha256'] != sha(plan['reference_plan'])
            or proof['status'] != 'passed_full_duplication_sister_reference_readback'
            or proof['producer_receipt_sha256'] != sha(reference / 'receipt.json')
            or Path(sourceplan['output']).resolve() != reference.resolve()):
        raise ValueError('Unverified native reference source')
    reviewplan = json.loads(Path(sourceplan['review_plan']).read_text())
    if Path(reviewplan['bridge']).resolve() != Path(sourceplan['bridge']).resolve():
        raise ValueError('Native model bridge differs')
    for name, digest in receipt['artifacts'].items():
        if sha(reference / name) != digest:
            raise ValueError('Changed native reference artifact')
    with sqlite3.connect('file:' + str(Path(sourceplan['bridge']).resolve()) + '?mode=ro', uri=True) as conn:
        models = {gene: (mid, version) for gene, mid, version in conn.execute(
            "SELECT taxon_id || '_' || protein_id,model_id,version FROM structures")}
    return reference, receipt, reviewplan, models


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pin: ' + path)

    verify()
    source, receipt, reviewplan, models = sources(plan)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    reports = []
    for entry in reviewplan['guides']:
        guide = entry['guide']
        groups = defaultdict(list)
        seen_rows = set()
        with (source / (guide + '_sister_references.tsv')).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = row['family'], row['gene_node'], row['gene_a'], row['gene_b']
                if key in seen_rows:
                    raise ValueError('Repeated candidate')
                seen_rows.add(key)
                groups[row['family']].append(row)
        counts, contexts, agreement, matrix, seen, writer = Counter(), Counter(), Counter(), Counter(), set(), None
        with (out / (guide + '_sequence_first_references.tsv')).open('w') as output, Path(entry['trees']).open() as handle:
            for line in handle:
                family, newick = line.rstrip().split(': ', 1)
                if family not in groups:
                    continue
                if family in seen:
                    raise ValueError('Repeated family tree')
                index = tree_index(newick)
                for original in groups[family]:
                    result = sequence_first(original, index, models)
                    counts[result['sequence_first_status']] += 1
                    contexts[original['status']] += 1
                    agreement[str(result['chosen_gene_matches_availability_choice'])] += 1
                    matrix[original['status'] + '|' + result['sequence_first_status']] += 1
                    data = {**original, **result}
                    for field in ['sequence_first_nearest_genes', 'sequence_first_tied_model_assignments']:
                        data[field] = json.dumps(data[field], separators=(',', ':'))
                    if writer is None:
                        writer = csv.DictWriter(output, list(data), delimiter='\t', lineterminator='\n')
                        writer.writeheader()
                    writer.writerow(data)
                seen.add(family)
                if len(seen) % 1000 == 0:
                    print(guide, 'sequence-first families', len(seen), flush=True)
        expected = next(r for r in receipt['guides'] if r['guide'] == guide)
        if seen != set(groups) or len(seen_rows) != expected['candidates'] or dict(contexts) != expected['counts']:
            raise ValueError('Incomplete source candidate universe')
        reports.append(dict(guide=guide, candidates=len(seen_rows), families=len(seen),
                            counts=dict(counts), source_context_counts=dict(contexts), chosen_gene_agreement=dict(agreement),
                            source_context_by_choice_status=dict(matrix)))
        print(json.dumps(reports[-1]), flush=True)
    verify()
    result = dict(status='complete_full_sequence_first_sister_reference_inventory_pending_readback',
                  plan_sha256=sha(args.plan), original_reference_receipt_sha256=sha(source / 'receipt.json'),
                  native_reference_readback_sha256=sha(plan['reference_readback']), guides=reports,
                  source_hashes=bindings, artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Every target retains original native availability-based reference fields. Nearest nonfocal sister genes selected solely by fixed native sequence-tree path before any model lookup; all ties within1e-12 and lexical chosen gene retained even when unmodeled. Availability-based distance increment and tie/model availability explicitly measured. Parent-ineligible targets retained, not promoted to valid references. Extant choices, not ancestors or independent orthology; no predictions, structural responses, alignment quality or biological asymmetry inference.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    csv.field_size_limit(32*1024*1024)
    main()
