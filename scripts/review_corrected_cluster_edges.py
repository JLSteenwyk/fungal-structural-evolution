#!/usr/bin/env python3
"""Apply verified inclusive-span score decisions to original representative edges."""
import argparse
import json
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['control_review', 'prior_review', 'output']:
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use immutable output')
    control = checked_receipt(a.control_review); prior = checked_receipt(a.prior_review)
    if control['original_completion_sha256'] != prior['exact_completion_sha256']:
        raise ValueError('Control and prior score sources differ')
    decisions = {}
    for row in read_table(a.control_review / 'directed_control_comparisons.tsv'):
        if row['stage'] != 'normalization':
            continue
        key = row['query'], row['target']
        if key in decisions or row['after_valid'] != 'True':
            raise ValueError('Duplicate pair or corrected score invalid')
        decisions[key] = row
    old = read_table(a.prior_review / 'edge_score_review.tsv')
    expected = {(r['representative_input_id'], r['member_input_id']) for r in old}
    if len(expected) != len(old) or set(decisions) != expected | {(b, a) for a, b in expected}:
        raise ValueError('Complete directed edge grid mismatch')
    rows = []; transitions = Counter()
    for row in old:
        rep, member = row['representative_input_id'], row['member_input_id']
        forward, reverse = decisions[rep, member], decisions[member, rep]
        old_passes = [forward['before_passes'] == 'True', reverse['before_passes'] == 'True']
        passes = [forward['after_passes'] == 'True', reverse['after_passes'] == 'True']
        status = lambda values: 'both_directions_pass' if all(values) else 'one_direction_passes' if any(values) else 'neither_direction_passes'
        if status(old_passes) != row['exact_status']:
            raise ValueError('Prior edge classification mismatch')
        new_status = status(passes)
        transitions[row['exact_status'] + '_to_' + new_status] += 1
        rows.append({'representative_input_id': rep, 'member_input_id': member,
                     'prior_exact_status': row['exact_status'], 'exact_status': new_status,
                     'forward_exact_passes': passes[0], 'reverse_exact_passes': passes[1],
                     'any_exact_score_out_of_bounds': False,
                     'normalization': 'inclusive_min_endpoint_span',
                     'membership_eligibility_changed': all(old_passes) != all(passes)})
    a.output.mkdir(parents=True); write_table(a.output / 'edge_score_review.tsv', rows)
    result = {'status': 'complete_corrected_normalization_bidirectional_edge_review',
              'validation_receipt_sha256': prior['validation_receipt_sha256'],
              'prior_review_receipt_sha256': sha(a.prior_review / 'receipt.json'),
              'control_review_receipt_sha256': sha(a.control_review / 'receipt.json'),
              'script_sha256': sha(Path(__file__)), 'directed_pairs': len(decisions), 'edge_pairs': len(rows),
              'edge_status_counts': dict(Counter(r['exact_status'] for r in rows)),
              'transitions': dict(transitions),
              'membership_eligibility_changes': sum(r['membership_eligibility_changed'] for r in rows),
              'interpretation': 'Fixed original representative edges under inclusive-span output normalization; existing exact_status fields retained for downstream schema compatibility. Both directions required. Preserves prior review and all source links. No new search, all-pairs equivalence, orthology or novelty claim.',
              'artifacts': {'edge_score_review.tsv': sha(a.output / 'edge_score_review.tsv')}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
