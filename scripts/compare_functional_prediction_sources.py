"""Compare full audited functional-site coverage, preserving predictor differences."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import pandas as pd

VARIABLE = ['source_label', 'model_id', 'model_version', 'model_available', 'site_status',
    'matrix_column_1based', 'paired_column_1based', 'focal_plddt', 'native_valid',
    'native_state', 'feature_min_plddt', 'feature_max_pae', 'joint_feature_confident',
    'observed_in_paired_alignment']
KEYS = ['marker', 'taxon_id', 'site_id']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    roots = dict(afdb=Path('results/functional_sites/afdb-recovered-linked-20260927-v1'),
                 esmfold=Path('results/functional_sites/esmfold-all-completed-linked-20260923-v1'))
    audits = dict(afdb=Path('metadata/recovered_afdb_functional_site_join_readback_20260927.json'),
                  esmfold=Path('metadata/completed_esmfold_functional_site_join_readback.json'))
    sources, tables, receipts = {}, {}, {}
    for label, root in roots.items():
        receipt = json.loads((root / 'receipt.json').read_text())
        audit = json.loads(audits[label].read_text())
        assert receipt['status'] == 'complete_functional_site_structure_alignment_join'
        assert audit['status'] == 'passed_full_completed_functional_site_row_readback'
        assert audit['producer_receipt_sha256'] == sha(root / 'receipt.json')
        path = root / 'site_structure_links.tsv'
        assert sha(path) == receipt['artifacts'][path.name]
        for p in (root / 'receipt.json', audits[label], path):
            sources[str(p)] = sha(p)
        rows = list(csv.DictReader(path.open(), delimiter='\t'))
        tables[label] = {tuple(row[k] for k in KEYS): row for row in rows}
        assert len(tables[label]) == len(rows) == 17105
        receipts[label] = receipt
    assert set(tables['afdb']) == set(tables['esmfold'])
    assert receipts['afdb']['source_receipts']['functional_sites'] == receipts['esmfold']['source_receipts']['functional_sites']
    shared = [k for k in next(iter(tables['afdb'].values())) if k not in VARIABLE]
    result = []
    for key in sorted(tables['afdb']):
        a, e = tables['afdb'][key], tables['esmfold'][key]
        assert set(a) == set(e)
        assert all(a[k] == e[k] for k in shared), key
        av, ev = a['observed_in_paired_alignment'] == 'True', e['observed_in_paired_alignment'] == 'True'
        for row, observed in ((a, av), (e, ev)):
            assert row['observed_in_paired_alignment'] in ('True', 'False')
            if observed:
                assert row['joint_feature_confident'] == 'True' and row['native_valid'] == 'True'
                assert row['protein_residue_1based'] and row['matrix_column_1based'] and row['paired_column_1based']
        if av and ev:
            assert a['matrix_column_1based'] == e['matrix_column_1based']
        coverage = 'both' if av and ev else 'afdb_only' if av else 'esmfold_only' if ev else 'neither'
        comparison = ('same' if a['native_state'] == e['native_state'] else 'different') if av and ev else 'not_comparable'
        row = {k: a[k] for k in shared}
        row.update({label + '_' + k: source[k] for label, source in [('afdb', a), ('esmfold', e)] for k in VARIABLE})
        row.update(paired_coverage=coverage, paired_state_comparison=comparison)
        result.append(row)
    # Aggregate repeated profile annotations only at exact protein coordinates.
    groups = defaultdict(list)
    residue_keys = ['marker', 'taxon_id', 'protein_id', 'protein_residue_1based']
    for row in result:
        if row['protein_residue_1based']:
            groups[tuple(row[k] for k in residue_keys)].append(row)
    residues = []
    consistent = ['observed_amino_acid', 'paired_coverage', 'paired_state_comparison']
    consistent += [p + '_' + k for p in roots for k in VARIABLE]
    for key, rows in sorted(groups.items()):
        assert all(len({r[k] for r in rows}) == 1 for k in consistent), key
        row = dict(zip(residue_keys, key))
        row.update({k: rows[0][k] for k in consistent})
        row.update(annotation_rows=len(rows),
                   any_conserved_candidate=any(r['conserved_candidate'] == 'True' for r in rows),
                   any_other_correspondence=any(r['conserved_candidate'] == 'False' for r in rows),
                   pfam_accessions=';'.join(sorted({r['pfam_accession'] for r in rows})))
        residues.append(row)
    args.output.mkdir(parents=True)
    write(args.output / 'annotation_source_comparison.tsv', result)
    write(args.output / 'residue_source_comparison.tsv', residues)
    # Read all exported source columns back through an independent dataframe join.
    frames = {label: pd.read_csv(root / 'site_structure_links.tsv', sep='\t', dtype=str, keep_default_na=False)
              for label, root in roots.items()}
    expected = frames['afdb'].merge(frames['esmfold'], on=shared, how='outer',
                                     suffixes=('_afdb', '_esmfold'), validate='one_to_one', indicator=True)
    assert expected['_merge'].eq('both').all() and len(expected) == 17105
    expected = expected.drop(columns='_merge').rename(columns={k+'_'+p: p+'_'+k for p in roots for k in VARIABLE})
    exported = pd.read_csv(args.output / 'annotation_source_comparison.tsv', sep='\t', dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(exported[expected.columns].sort_values(KEYS).reset_index(drop=True),
                                  expected.sort_values(KEYS).reset_index(drop=True))
    a = exported.afdb_observed_in_paired_alignment.eq('True')
    e = exported.esmfold_observed_in_paired_alignment.eq('True')
    masks = {'both': a & e, 'afdb_only': a & ~e, 'esmfold_only': ~a & e, 'neither': ~a & ~e}
    for label, mask in masks.items():
        assert exported.paired_coverage.eq(label).equals(mask)
    different = a & e & exported.afdb_native_state.ne(exported.esmfold_native_state)
    assert exported.paired_state_comparison.eq('different').equals(different)
    assert exported.paired_state_comparison.eq('same').equals(a & e & ~different)
    rd = pd.read_csv(args.output / 'residue_source_comparison.tsv', sep='\t', dtype=str, keep_default_na=False)
    valid = exported[exported.protein_residue_1based.ne('')]
    grouping = valid.groupby(residue_keys)
    expected_rd = grouping[consistent].first()
    expected_rd['annotation_rows'] = grouping.size().astype(str)
    expected_rd['any_conserved_candidate'] = grouping.conserved_candidate.apply(lambda x: str(x.eq('True').any()))
    expected_rd['any_other_correspondence'] = grouping.conserved_candidate.apply(lambda x: str(x.eq('False').any()))
    expected_rd['pfam_accessions'] = grouping.pfam_accession.apply(lambda x: ';'.join(sorted(set(x))))
    pd.testing.assert_frame_equal(rd.set_index(residue_keys).sort_index(), expected_rd[rd.columns.difference(residue_keys, sort=False)].sort_index())
    for path, digest in sources.items():
        assert sha(path) == digest
    report = dict(status='complete_full_functional_prediction_source_comparison_with_readback',
        annotation_rows=len(result), coordinate_residues=len(residues),
        annotation_coverage_counts=dict(Counter(r['paired_coverage'] for r in result)),
        residue_coverage_counts=dict(Counter(r['paired_coverage'] for r in residues)),
        annotation_state_comparison_counts=dict(Counter(r['paired_state_comparison'] for r in result)),
        residue_state_comparison_counts=dict(Counter(r['paired_state_comparison'] for r in residues)),
        source_sha256=sources, script_sha256=sha(__file__),
        artifacts={p.name: sha(p) for p in args.output.iterdir()},
        scope='Complete annotation universe and exact-coordinate deduplication. Both predictor columns preserved and exported rows independently joined/read back. Conditional observed-state agreement is not accuracy, independent validation, evolutionary change, catalytic function or a calibrated predictor-error rate. Structural-alphabet states depend on nonlocal context; missing coverage is not a state difference.')
    (args.output / 'receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('source_sha256', 'artifacts')}, indent=2))


if __name__ == '__main__':
    main()
