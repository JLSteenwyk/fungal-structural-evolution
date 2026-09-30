#!/usr/bin/env python3
"""Independently outer-merge all sequence-first guide choices and check every export."""
import argparse
import json
from pathlib import Path

import pandas as pd
from screen_duplication_alignment_reuse import sha


def load(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def expected(left, right):
    fields = ['family', 'status', 'sequence_first_status', 'sequence_first_nearest_genes',
              'sequence_first_chosen_gene', 'sequence_first_chosen_model', 'sequence_first_chosen_version',
              'sequence_first_minimum_distance', 'availability_distance_increment']
    frames = []
    for guide, original in [('profile', left), ('mafft', right)]:
        data = original.copy()
        keys = [sorted(pair) for pair in zip(data.gene_a, data.gene_b)]
        data['gene_a'] = [k[0] for k in keys]
        data['gene_b'] = [k[1] for k in keys]
        if data.duplicated(['gene_a', 'gene_b']).any():
            raise ValueError('Duplicate input pair')
        for text, chosen in zip(data.sequence_first_nearest_genes, data.sequence_first_chosen_gene):
            ties = json.loads(text)
            if ties != sorted(set(ties)) or chosen != (ties[0] if ties else ''):
                raise ValueError('Invalid source tie/choice')
        frames.append(data[['gene_a', 'gene_b', *fields]].rename(columns={f: guide + '_' + f for f in fields}))
    out = frames[0].merge(frames[1], how='outer', on=['gene_a', 'gene_b'], validate='one_to_one', indicator=True)
    out['_merge'] = out['_merge'].astype(str)
    shared = out['_merge'] == 'both'
    out = out.fillna('')
    for guide, side in [('profile', 'right_only'), ('mafft', 'left_only')]:
        absent = out['_merge'] == side
        out.loc[absent, [guide + '_status', guide + '_sequence_first_status']] = 'not_candidate'
    eligible = shared.copy()
    for guide in ['profile', 'mafft']:
        eligible &= out[guide + '_status'].isin(['no_modeled_nonfocal_sister', 'provisional_reference_available'])
    relations = []
    for p, m, kind, ok in zip(out.profile_sequence_first_nearest_genes,
                              out.mafft_sequence_first_nearest_genes, out['_merge'], eligible):
        if kind == 'left_only':
            value = 'profile_only_candidate'
        elif kind == 'right_only':
            value = 'mafft_only_candidate'
        elif not ok:
            value = 'ineligible_sequence_context_in_one_or_both'
        else:
            a, b = frozenset(json.loads(p)), frozenset(json.loads(m))
            if not a or not b:
                raise ValueError('Eligible source lacks sister genes')
            value = ('same_nearest_sequence_gene_set' if a == b else
                     'disjoint_nearest_sequence_gene_sets' if a.isdisjoint(b) else
                     'overlapping_nearest_sequence_gene_sets')
        relations.append(value)
    out['relation'] = relations
    same_gene = out.profile_sequence_first_chosen_gene == out.mafft_sequence_first_chosen_gene
    model_same = ((out.profile_sequence_first_chosen_model == out.mafft_sequence_first_chosen_model) &
                  (out.profile_sequence_first_chosen_version == out.mafft_sequence_first_chosen_version))
    available = eligible & (out.profile_sequence_first_chosen_model != '') & (out.mafft_sequence_first_chosen_model != '')
    if (eligible & same_gene & ~model_same).any():
        raise ValueError('Same gene has conflicting frozen model/version')
    out['same_lexical_gene'] = same_gene.astype(int).astype(str).where(eligible, '')
    out['both_lexical_models_available'] = available.astype(int).astype(str).where(eligible, '')
    out['same_lexical_model'] = model_same.astype(int).astype(str).where(available, '')
    totals = dict(union_pairs=len(out), shared_pairs=int(shared.sum()),
                  sequence_eligible_in_both=int(eligible.sum()), same_lexical_gene=int((eligible & same_gene).sum()),
                  both_lexical_models_available=int(available.sum()), same_lexical_model=int((available & model_same).sum()),
                  relations=out.relation.value_counts().to_dict())
    matrices = {}
    for name, field in [('context_matrix.tsv', 'status'), ('choice_matrix.tsv', 'sequence_first_status')]:
        matrix = out.groupby(['profile_' + field, 'mafft_' + field], sort=True).size().reset_index(name='pairs')
        matrix['pairs'] = matrix.pairs.astype(str)
        matrices[name] = matrix
    return out.drop(columns='_merge').sort_values(['gene_a', 'gene_b']).reset_index(drop=True), matrices, totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    root, source = Path(plan['output']), Path(plan['inventory'])
    rp, nrp = root / 'receipt.json', source / 'receipt.json'
    record, native = json.loads(rp.read_text()), json.loads(nrp.read_text())
    proof = json.loads(Path(plan['native_readback']).read_text())
    if (record['status'] != 'complete_full_sequence_first_guide_comparison_pending_readback'
            or record['plan_sha256'] != sha(args.plan)
            or native['status'] != 'complete_full_sequence_first_sister_reference_inventory_pending_readback'
            or native['plan_sha256'] != sha(plan['source_plan'])
            or proof['status'] != 'passed_full_sequence_first_sister_reference_native_readback'
            or proof['producer_receipt_sha256'] != sha(nrp)
            or record['native_receipt_sha256'] != sha(nrp)
            or record['native_readback_sha256'] != sha(plan['native_readback'])):
        raise ValueError('Unverified comparison/native source')
    bindings = {str(args.plan): sha(args.plan), **plan['pins'], str(rp): sha(rp)}
    for folder, receipt in [(root, record), (source, native)]:
        for name, digest in receipt['artifacts'].items():
            path = str(folder / name)
            if path in bindings and bindings[path] != digest:
                raise ValueError('Conflicting artifact binding')
            bindings[path] = digest
    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed pinned source/export: ' + path)
    verify()
    frames = [load(source / (g + '_sequence_first_references.tsv')) for g in ['profile', 'mafft']]
    for guide, frame in zip(['profile', 'mafft'], frames):
        if len(frame) != next(r['candidates'] for r in native['guides'] if r['guide'] == guide):
            raise ValueError('Incomplete source universe')
    table, matrices, totals = expected(*frames)
    pd.testing.assert_frame_equal(load(root / 'sequence_first_guide_comparison.tsv'), table, check_like=False)
    for name, matrix in matrices.items():
        pd.testing.assert_frame_equal(load(root / name), matrix, check_like=False)
    if any(record[k] != value for k, value in totals.items()):
        raise ValueError('Exported totals differ')
    verify()
    result = dict(status='passed_full_sequence_first_guide_comparison_readback',
                  producer_receipt_sha256=sha(rp), source_hashes=bindings, **totals,
                  checker_sha256=sha(__file__), scope='Full independent dataframe outer merge checks every '
                  'gene-pair choice, tie, source context, missing assignment, matrix and summary. Requires '
                  'full native sequence-first choice proof; no structure availability substitution, '
                  'orthology acceptance, structural geometry or biological asymmetry inference.')
    with args.output.open('x') as handle:
        handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', *totals]}), flush=True)


if __name__ == '__main__':
    main()
