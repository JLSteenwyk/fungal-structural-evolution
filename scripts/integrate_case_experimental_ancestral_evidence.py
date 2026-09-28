#!/usr/bin/env python3
"""Preserve all experimental sensitivity strata beside all ancestral case evidence."""
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    pins = {}

    def read(root, name):
        root = Path(root)
        receipt = root / 'receipt.json'
        data = root / name
        assert sha(data) == json.loads(receipt.read_text())['artifacts'][name]
        pins.update({str(p): sha(p) for p in [receipt, data]})
        with data.open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    cases = read('results/structural_comparisons/case-ancestral-uncertainty-integration-20260927-v1',
                 'case_dossiers_with_ancestral_uncertainty.tsv')
    root = Path('results/experimental_structures/whole-domain-case-reference-robustness-20260927-v1')
    closure = Path('metadata/case_experimental_reference_robustness_completed_20260927.json')
    c = json.loads(closure.read_text())
    assert sha(root / 'receipt.json') == c['producer_receipt_sha256']
    assert sha(c['readback_receipt_path']) == c['readback_receipt_sha256']
    pins[str(closure)] = sha(closure)
    strata = read(root, 'case_summary.tsv')
    keys = ['family', 'gene_a', 'gene_b', 'pfam_accession', 'species_name']
    key = lambda row: tuple(row[k] for k in keys)
    index = {key(row): row for row in cases}
    assert len(index) == len(cases) == 13 and len(strata) == 936
    assert set(index) == {key(row) for row in strata}
    axes = ['metric', 'screen', 'margin_angstrom']
    expected = {(r['metric'], r['screen'], r['margin_angstrom']) for r in strata}
    assert len(expected) == 72
    for k in index:
        rows = [r for r in strata if key(r) == k]
        assert len(rows) == 72
        assert {tuple(r[a] for a in axes) for r in rows} == expected
    output = []
    for row in strata:
        assert int(row['qualified_units']) == sum(int(row[x]) for x in
            ['same_direction_units', 'opposite_direction_units', 'variable_or_within_margin_units'])
        assert 0 <= int(row['qualified_units']) <= int(row['candidate_units'])
        joined = dict(index[key(row)])
        joined.update({'experimental_' + k: v for k, v in row.items() if k not in keys})
        output.append(joined)
    out = Path('results/structural_comparisons/case-experimental-ancestral-evidence-20260927-v1')
    out.mkdir(exist_ok=False)
    target = out / 'case_evidence_all_settings.tsv'
    with target.open('w') as stream:
        writer = csv.DictWriter(stream, list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    # A separate relational merge checks every emitted value and all membership.
    experimental = pd.DataFrame(strata).rename(columns={k: 'experimental_' + k for k in strata[0] if k not in keys})
    expected_frame = experimental.merge(pd.DataFrame(cases), on=keys, how='left', validate='many_to_one')
    actual = pd.read_csv(target, sep='\t', dtype=str, keep_default_na=False)
    sort_keys = keys + ['experimental_' + a for a in axes]
    pd.testing.assert_frame_equal(actual.sort_values(sort_keys).reset_index(drop=True),
        expected_frame[list(actual.columns)].sort_values(sort_keys).reset_index(drop=True))
    receipt = dict(status='complete_all_case_experimental_ancestral_integration',
        cases=13, settings_per_case=72, rows=len(output), fields=len(output[0]),
        zero_qualified_rows=sum(int(r['qualified_units']) == 0 for r in strata),
        pins=pins, script_sha256=sha(__file__), artifacts={target.name: sha(target)},
        verification='Every emitted field checked against an independent pandas many-to-one merge; all 13 by 72 strata preserved.',
        scope='Descriptive integration only. Reference units are dependent, ancestral counts repeat fits and nodes. No pooling across settings, case ranking, significance, mechanism, ancestral polarity or qualified ancestral ensemble is implied.')
    rp = out / 'receipt.json'
    rp.write_text(json.dumps(receipt, indent=2) + '\n')
    receipt.update(completed_receipt_path=str(rp), completed_receipt_sha256=sha(rp))
    Path('metadata/case_experimental_ancestral_integration_completed_20260927.json').write_text(json.dumps(receipt, indent=2) + '\n')
    Path('docs/tables/case_evidence_all_settings_20260927.tsv').write_bytes(target.read_bytes())
    print(json.dumps({k: receipt[k] for k in ['cases', 'rows', 'fields', 'zero_qualified_rows']}))


if __name__ == '__main__':
    main()
