#!/usr/bin/env python3
"""Add complete local tree/fit diagnostics while retaining historical warnings."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    with Path(p).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == len({r['case_id'] for r in rows})
    return {r['case_id']: r for r in rows}


def main():
    base = Path('results/cds/full-group-alignment-readiness-20260927-v1')
    tree = Path('results/phylogeny/local-codon-tree-audit-20260927-v1')
    fit = Path('results/cds/local-mg94-audit-20260927-v1')
    inputs = [base/'cases.tsv', base/'receipt.json', tree/'receipt.json', fit/'receipt.json', Path('metadata/local_mg94_completed_audit_20260927.json'), Path('metadata/local_codon_tree_information_20260927.tsv')]
    for root in (base, tree, fit):
        r = json.loads((root/'receipt.json').read_text())
        for name, digest in r['artifacts'].items():
            assert sha(root/name) == digest
    tproof = json.loads((tree/'receipt.json').read_text())
    fproof = json.loads((fit/'receipt.json').read_text())
    assert tproof['status'] == 'passed_full_genus_tree_audit'
    assert fproof['status'] == 'passed_saved_fit_readback'
    completion = json.loads(inputs[4].read_text())
    assert completion['audit_receipt_sha256'] == sha(fit/'receipt.json')
    pins = {str(p): sha(p) for p in inputs}
    ledger, trees, fits, info = read(base/'cases.tsv'), read(tree/'case_audit.tsv'), read(fit/'case_audit.tsv'), read(inputs[-1])
    assert set(ledger) == set(info) and len(ledger) == 1712
    assert set(trees) == set(fits) and len(fits) == 1632
    out = []
    counts = Counter()
    for case, old in ledger.items():
        r = dict(old)
        r['historical_flags_scope'] = 'original_alignment_diagnostics_preserved;local_diagnostics_in_separate_columns'
        r['local_fit_completed'] = case in fits
        r['local_tree_disposition'] = info[case]['status']
        flags = []
        if case in fits:
            t, f = trees[case], fits[case]
            assert t['taxa'] == f['taxa'] == info[case]['fcs_remaining_taxa']
            assert int(t['nucleotide_columns']) == 3*int(f['codons']) == int(old['local_nucleotide_columns'])
            assert t['marker_copy_caveat'] == f['marker_copy_caveat'] == old['copy_caveat']
            for column in ('near_zero_edges_le1e_5', 'long_edges_ge10', 'warning_lines', 'gamma_alpha', 'total_tree_length'):
                r['local_tree_'+column] = t[column]
            r['local_fit_numerical_review_flags'] = f['numerical_review_flags']
            r['local_fit_warning_lines'] = f['warning_lines']
            if int(t['near_zero_edges_le1e_5']): flags.append('near_zero_nucleotide_tree_edge')
            if int(t['long_edges_ge10']): flags.append('very_long_nucleotide_tree_edge')
            if int(t['warning_lines']): flags.append('tree_warning_review_required')
            if f['numerical_review_flags']: flags.append('saved_fit_numerical_review')
            if int(f['warning_lines']): flags.append('fit_warning_review_required')
        else:
            for column in ('near_zero_edges_le1e_5', 'long_edges_ge10', 'warning_lines', 'gamma_alpha', 'total_tree_length'):
                r['local_tree_'+column] = ''
            r['local_fit_numerical_review_flags'] = ''
            r['local_fit_warning_lines'] = ''
            flags.append('local_fit_unavailable')
        if old['copy_caveat'] != 'none_recorded': flags.append('copy_reconciliation_required')
        r['local_case_review_flags'] = ';'.join(flags)
        r['selection_eligibility'] = 'not_established'
        counts.update(flags)
        out.append(r)
    root = Path('results/cds/local-codon-diagnostic-review-20260927-v1')
    root.mkdir(parents=True, exist_ok=False)
    table = root/'cases.tsv'
    with table.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(out)
    for path, digest in pins.items(): assert sha(path) == digest
    result = dict(status='complete_local_codon_diagnostic_review_pending_readback', cases=len(out), fitted_cases=len(fits), local_flag_counts=dict(counts), historical_flagged_cases=sum(bool(r['case_specific_review_flags']) for r in out), local_flagged_cases=sum(bool(r['local_case_review_flags']) for r in out), selection_eligible_cases=0, sources=pins, script_sha256=sha(__file__), artifacts={'cases.tsv':sha(table)}, scope='Full local diagnostics joined to unchanged historical warnings. Existing near-zero (<=1e-5) and long-edge (>=10) labels are review triggers, not biological exclusion thresholds. Warning-line counts are not numbers of independent problems. Missing fits remain missing. No historical concern is cleared and no selection eligibility established.')
    (root/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
