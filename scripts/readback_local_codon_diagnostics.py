#!/usr/bin/env python3
"""Reconstruct the full diagnostic join with independent indexed dataframes."""
import hashlib
import json
from pathlib import Path
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    frame = pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False).set_index('case_id', verify_integrity=True)
    return frame.sort_index()


def main():
    root = Path('results/cds/local-codon-diagnostic-review-20260927-v1')
    receipt = json.loads((root/'receipt.json').read_text())
    for p, digest in receipt['sources'].items(): assert sha(p) == digest
    assert sha(root/'cases.tsv') == receipt['artifacts']['cases.tsv']
    actual = read(root/'cases.tsv')
    old = read('results/cds/full-group-alignment-readiness-20260927-v1/cases.tsv')
    tree = read('results/phylogeny/local-codon-tree-audit-20260927-v1/case_audit.tsv')
    fit = read('results/cds/local-mg94-audit-20260927-v1/case_audit.tsv')
    info = read('metadata/local_codon_tree_information_20260927.tsv')
    assert tree.index.equals(fit.index) and old.index.equals(actual.index) and old.index.equals(info.index)
    columns = old.columns.difference(['historical_flags_scope'])
    pd.testing.assert_frame_equal(actual[columns], old[columns])
    assert (actual.historical_flags_scope == 'original_alignment_diagnostics_preserved;local_diagnostics_in_separate_columns').all()
    assert actual.local_fit_completed.tolist() == ['True' if case in fit.index else 'False' for case in old.index]
    assert actual.local_tree_disposition.equals(info.status.rename('local_tree_disposition'))
    for col in ('near_zero_edges_le1e_5','long_edges_ge10','warning_lines','gamma_alpha','total_tree_length'):
        assert actual['local_tree_'+col].tolist() == tree[col].reindex(old.index, fill_value='').tolist()
    for col in ('numerical_review_flags','warning_lines'):
        assert actual['local_fit_'+col].tolist() == fit[col].reindex(old.index, fill_value='').tolist()
    flags = pd.DataFrame(index=old.index)
    flags['near_zero_nucleotide_tree_edge'] = tree.near_zero_edges_le1e_5.astype(int).gt(0).reindex(old.index,fill_value=False)
    flags['very_long_nucleotide_tree_edge'] = tree.long_edges_ge10.astype(int).gt(0).reindex(old.index,fill_value=False)
    flags['tree_warning_review_required'] = tree.warning_lines.astype(int).gt(0).reindex(old.index,fill_value=False)
    flags['saved_fit_numerical_review'] = fit.numerical_review_flags.ne('').reindex(old.index,fill_value=False)
    flags['fit_warning_review_required'] = fit.warning_lines.astype(int).gt(0).reindex(old.index,fill_value=False)
    flags['local_fit_unavailable'] = ~old.index.isin(fit.index)
    flags['copy_reconciliation_required'] = old.copy_caveat.ne('none_recorded')
    for case, r in flags.iterrows():
        assert set(actual.loc[case,'local_case_review_flags'].split(';')) - {''} == set(r.index[r])
    counts = {k:int(v) for k,v in flags.sum().items() if v}
    assert counts == receipt['local_flag_counts']
    assert receipt['local_flagged_cases'] == int(flags.any(axis=1).sum())
    assert receipt['historical_flagged_cases'] == int(old.case_specific_review_flags.ne('').sum())
    assert receipt['cases'] == len(old) == 1712 and receipt['fitted_cases'] == len(fit) == 1632
    assert (actual.selection_eligibility == 'not_established').all() and receipt['selection_eligible_cases'] == 0
    proof = dict(status='passed_full_local_codon_diagnostic_join_readback',cases=len(old),fitted_cases=len(fit),local_flag_counts=counts,source_receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(__file__),scope='Independent dataframe reconstruction verifies every original field, local diagnostic column, disposition, review flag and aggregate. Does not establish biological adequacy or clear historical concerns.')
    Path('metadata/local_codon_diagnostic_review_readback_20260927.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))


if __name__ == '__main__':
    main()
