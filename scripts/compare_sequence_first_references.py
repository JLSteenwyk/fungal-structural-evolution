#!/usr/bin/env python3
"""Compare the full sequence-first sister-reference universe across native guides."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from screen_duplication_alignment_reuse import sha


ELIGIBLE = {'provisional_reference_available', 'no_modeled_nonfocal_sister'}
SOURCE_FIELDS = ['family', 'status', 'sequence_first_status',
                 'sequence_first_nearest_genes', 'sequence_first_chosen_gene',
                 'sequence_first_chosen_model', 'sequence_first_chosen_version',
                 'sequence_first_minimum_distance', 'availability_distance_increment']
FIELDS = ['gene_a', 'gene_b'] + [g + '_' + f for g in ['profile', 'mafft'] for f in SOURCE_FIELDS] + [
    'relation', 'same_lexical_gene', 'both_lexical_models_available', 'same_lexical_model']


def source_bindings(plan):
    root = Path(plan['inventory'])
    rp = root / 'receipt.json'
    receipt = json.loads(rp.read_text())
    proof = json.loads(Path(plan['native_readback']).read_text())
    producer_plan = json.loads(Path(plan['source_plan']).read_text())
    if (receipt['status'] != 'complete_full_sequence_first_sister_reference_inventory_pending_readback'
            or receipt['plan_sha256'] != sha(plan['source_plan'])
            or Path(producer_plan['output']).resolve() != root.resolve()
            or proof['status'] != 'passed_full_sequence_first_sister_reference_native_readback'
            or proof['producer_receipt_sha256'] != sha(rp)):
        raise ValueError('Unverified sequence-first source')
    bindings = dict(plan['pins'])
    for path, digest in {str(rp): sha(rp), str(plan['native_readback']): sha(plan['native_readback']),
                         str(plan['source_plan']): sha(plan['source_plan']),
                         **{str(root / name): digest for name, digest in receipt['artifacts'].items()}}.items():
        if path not in bindings or bindings[path] != digest:
            raise ValueError('Missing/conflicting source pin: ' + path)
    return receipt, bindings


def validate_row(row):
    ties = json.loads(row['sequence_first_nearest_genes'])
    if ties != sorted(set(ties)) or row['sequence_first_chosen_gene'] != (ties[0] if ties else ''):
        raise ValueError('Invalid sequence-first tie/lexical representation')
    if bool(row['sequence_first_chosen_model']) != bool(row['sequence_first_chosen_version']):
        raise ValueError('Partial model/version assignment')
    if row['status'] in ELIGIBLE and not ties:
        raise ValueError('Eligible sequence context lacks nonfocal sister genes')


def comparison(key, left, right):
    output = dict(zip(['gene_a', 'gene_b'], key))
    for guide, row in [('profile', left), ('mafft', right)]:
        output.update({guide + '_' + field: row[field] if row else
                       ('not_candidate' if field in ['status', 'sequence_first_status'] else '')
                       for field in SOURCE_FIELDS})
    if left is None:
        relation = 'mafft_only_candidate'
    elif right is None:
        relation = 'profile_only_candidate'
    elif left['status'] not in ELIGIBLE or right['status'] not in ELIGIBLE:
        relation = 'ineligible_sequence_context_in_one_or_both'
    else:
        a, b = [set(json.loads(row['sequence_first_nearest_genes'])) for row in [left, right]]
        relation = ('same_nearest_sequence_gene_set' if a == b else
                    'overlapping_nearest_sequence_gene_sets' if a & b else
                    'disjoint_nearest_sequence_gene_sets')
    output['relation'] = relation
    eligible = relation in {'same_nearest_sequence_gene_set', 'overlapping_nearest_sequence_gene_sets',
                            'disjoint_nearest_sequence_gene_sets'}
    output.update(same_lexical_gene='', both_lexical_models_available='', same_lexical_model='')
    if eligible:
        same = left['sequence_first_chosen_gene'] == right['sequence_first_chosen_gene']
        available = all(row['sequence_first_chosen_model'] for row in [left, right])
        output['same_lexical_gene'] = str(int(same))
        output['both_lexical_models_available'] = str(int(available))
        if available:
            output['same_lexical_model'] = str(int(all(left[field] == right[field] for field in
                                                     ['sequence_first_chosen_model', 'sequence_first_chosen_version'])))
        if same and any(left[field] != right[field] for field in
                        ['sequence_first_chosen_model', 'sequence_first_chosen_version']):
            raise ValueError('Same lexical gene has conflicting frozen model/version')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    receipt, bindings = source_bindings(plan)
    bindings[str(args.plan)] = sha(args.plan)
    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)
    verify()
    tables = {}
    for guide in ['profile', 'mafft']:
        rows = {}
        with (Path(plan['inventory']) / (guide + '_sequence_first_references.tsv')).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                validate_row(row)
                key = tuple(sorted([row['gene_a'], row['gene_b']]))
                if key in rows or key[0] == key[1]:
                    raise ValueError('Repeated/invalid duplicate gene pair')
                rows[key] = row
        if len(rows) != next(r['candidates'] for r in receipt['guides'] if r['guide'] == guide):
            raise ValueError('Source candidate scope differs')
        tables[guide] = rows
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    relations, status_matrix, choice_matrix = Counter(), Counter(), Counter()
    totals = Counter()
    with (out / 'sequence_first_guide_comparison.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key in sorted(set(tables['profile']) | set(tables['mafft'])):
            row = comparison(key, tables['profile'].get(key), tables['mafft'].get(key))
            writer.writerow(row)
            relations[row['relation']] += 1
            status_matrix[row['profile_status'], row['mafft_status']] += 1
            choice_matrix[row['profile_sequence_first_status'], row['mafft_sequence_first_status']] += 1
            totals['union_pairs'] += 1
            totals['shared_pairs'] += int(key in tables['profile'] and key in tables['mafft'])
            totals['sequence_eligible_in_both'] += int(row['same_lexical_gene'] != '')
            totals['same_lexical_gene'] += int(row['same_lexical_gene'] == '1')
            totals['both_lexical_models_available'] += int(row['both_lexical_models_available'] == '1')
            totals['same_lexical_model'] += int(row['same_lexical_model'] == '1')
    for name, matrix, suffix in [('context_matrix.tsv', status_matrix, 'status'),
                                 ('choice_matrix.tsv', choice_matrix, 'sequence_first_status')]:
        with (out / name).open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['profile_' + suffix, 'mafft_' + suffix, 'pairs'])
            writer.writerows((*key, value) for key, value in sorted(matrix.items()))
    verify()
    result = dict(status='complete_full_sequence_first_guide_comparison_pending_readback',
                  plan_sha256=sha(args.plan), native_receipt_sha256=sha(Path(plan['inventory']) / 'receipt.json'),
                  native_readback_sha256=sha(plan['native_readback']), **totals, relations=dict(relations),
                  source_hashes=bindings, artifacts={p.name: sha(p) for p in out.iterdir()},
                  scope='Full unordered gene-pair union, including every missing model, tie, parent-ineligible '
                  'context and one-guide candidate. Sequence eligibility fixed by original parent context, '
                  'not model availability. Empty model assignments never count as identical predicted models. '
                  'Guide sensitivity only; no accepted orthology, topology, ancestry or structural effect.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', *totals, 'relations']}), flush=True)


if __name__ == '__main__':
    csv.field_size_limit(32 * 1024 * 1024)
    main()
