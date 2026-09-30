#!/usr/bin/env python3
"""Independently join every original reference choice across both guide methods."""
import argparse
import json
from pathlib import Path

import pandas as pd

from run_ortholog_pair_guide_comparison import sha


def load(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def keyed(frame):
    frame = frame.copy()
    pairs = [tuple(sorted(pair)) for pair in zip(frame['gene_a'], frame['gene_b'])]
    frame['a'] = [pair[0] for pair in pairs]
    frame['b'] = [pair[1] for pair in pairs]
    if frame.duplicated(['a', 'b']).any():
        raise ValueError('Repeated source gene pair')
    return frame.set_index(['a', 'b']).sort_index()


def compare_sources(left, right):
    left, right = keyed(left), keyed(right)
    fields = ['status', 'chosen_reference_gene', 'reference_model', 'reference_version', 'nearest_reference_genes']
    union = left[fields].join(right[fields], how='outer', lsuffix='_profile', rsuffix='_mafft')
    present_left = union.index.isin(left.index)
    present_right = union.index.isin(right.index)
    both = present_left & present_right
    union = union.fillna('')
    output = pd.DataFrame(index=union.index)
    output['profile_status'] = union['status_profile'].where(present_left, 'not_candidate')
    output['mafft_status'] = union['status_mafft'].where(present_right, 'not_candidate')
    eligible = (output['profile_status'] == 'provisional_reference_available') & (output['mafft_status'] == 'provisional_reference_available')
    output['relation'] = 'not_provisional_in_both'
    output.loc[~present_right, 'relation'] = 'profile_only_candidate'
    output.loc[~present_left, 'relation'] = 'mafft_only_candidate'
    for key in union.index[eligible]:
        a = frozenset(json.loads(union.loc[key, 'nearest_reference_genes_profile']))
        b = frozenset(json.loads(union.loc[key, 'nearest_reference_genes_mafft']))
        if not a or not b:
            raise ValueError('Eligible reference set is empty')
        relation = ('same_nearest_reference_set' if a == b else
                    'disjoint_nearest_reference_sets' if a.isdisjoint(b) else 'overlapping_nearest_reference_sets')
        output.loc[key, 'relation'] = relation
    output['profile_reference'] = union['chosen_reference_gene_profile']
    output['mafft_reference'] = union['chosen_reference_gene_mafft']
    same_choice = union['chosen_reference_gene_profile'] == union['chosen_reference_gene_mafft']
    same_model = ((union['reference_model_profile'] == union['reference_model_mafft']) &
                  (union['reference_version_profile'] == union['reference_version_mafft']))
    if (eligible & same_choice & ~same_model).any():
        raise ValueError('Same chosen reference gene has different frozen model/version')
    output['same_chosen_reference'] = same_choice.astype(int).astype(str).where(eligible, '')
    output['same_chosen_model'] = same_model.astype(int).astype(str).where(eligible, '')
    output = output.reset_index().rename(columns={'a': 'gene_a', 'b': 'gene_b'})
    summary = dict(shared_duplicate_pairs=int(both.sum()), provisional_in_both=int(eligible.sum()),
                   same_chosen_reference_gene=int((eligible & same_choice).sum()),
                   same_chosen_reference_model=int((eligible & same_model).sum()),
                   relations=output['relation'].value_counts().to_dict())
    matrix = output.groupby(['profile_status', 'mafft_status'], sort=True).size().reset_index(name='pairs')
    matrix['pairs'] = matrix['pairs'].astype(str)
    return output, matrix, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'inventory', 'reference-readback', 'review', 'output']:
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    receipt_path = args.source / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    native_path = args.inventory / 'receipt.json'
    native = json.loads(native_path.read_text())
    audit = json.loads(args.reference_readback.read_text())
    review_path = args.review / 'receipt.json'
    review = json.loads(review_path.read_text())
    if (receipt['status'] != 'complete_duplication_reference_guide_comparison'
            or native['status'] != 'complete_duplication_sister_reference_inventory'
            or audit['status'] != 'passed_full_duplication_sister_reference_readback'
            or audit['producer_receipt_sha256'] != sha(native_path)):
        raise ValueError('Unverified guide/reference source')
    bindings = {**receipt['source_pins'], str(receipt_path): sha(receipt_path),
                str(args.reference_readback): sha(args.reference_readback)}
    for path in [native_path, review_path]:
        if bindings.get(str(path)) != sha(path):
            raise ValueError('Changed CLI source: ' + str(path))
    for root, record in [(args.source, receipt), (args.inventory, native),
                         (args.review, review), (args.reference_readback.parent, audit)]:
        for name, digest in record.get('artifacts', {}).items():
            path = str(root / name)
            if path in bindings and bindings[path] != digest:
                raise ValueError('Conflicting artifact binding')
            bindings[path] = digest

    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pinned input/output: ' + path)

    verify()
    frames = []
    source_counts = {}
    for guide in ['profile', 'mafft']:
        frame = load(args.inventory / (guide + '_sister_references.tsv'))
        original = keyed(load(args.review / (guide + '_candidate_tree_checks.tsv')))
        actual = keyed(frame)
        pd.testing.assert_frame_equal(actual[original.columns], original, check_like=False)
        source_counts[guide] = len(frame)
        if len(frame) != next(r for r in native['guides'] if r['guide'] == guide)['candidates']:
            raise ValueError('Full source candidate scope differs')
        frames.append(frame)
    expected, matrix, summary = compare_sources(*frames)
    pd.testing.assert_frame_equal(load(args.source / 'pair_reference_comparison.tsv'), expected, check_like=False)
    pd.testing.assert_frame_equal(load(args.source / 'eligibility_matrix.tsv'), matrix, check_like=False)
    if any(receipt[name] != value for name, value in summary.items()):
        raise ValueError('Guide summary differs')
    if summary['shared_duplicate_pairs'] != review['overlap']['shared_pairs']:
        raise ValueError('Candidate overlap differs')
    verify()
    result = dict(status='passed_full_reference_guide_comparison_readback',
                  producer_receipt_sha256=sha(receipt_path), reference_readback_sha256=sha(args.reference_readback),
                  source_guides=source_counts, union_pairs=len(expected), **summary,
                  source_hashes=bindings, checker_sha256=sha(__file__),
                  scope='Independent full dataframe outer join reproduces every saved reference-choice and '
                  'eligibility-matrix row, original candidate field and summary across both guides. All ties, '
                  'one-guide targets and unavailable references retained. Requires separate full native-tree '
                  'choice/path readback; not independent biological orthology, ancestral states, topology '
                  'acceptance or structural asymmetry inference.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
