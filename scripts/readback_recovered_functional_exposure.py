"""Reconstruct the complete functional/exposure join without producer helpers."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import pandas as pd


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    opts = parser.parse_args()
    plan = json.loads(opts.plan.read_text())
    ph = sha(opts.plan)
    for path, expected in plan['pins'].items():
        assert sha(path) == expected, path
    args = {key: Path(value) for key, value in plan['arguments'].items()}
    root = args['output']
    receipt = json.loads((root / 'receipt.json').read_text())
    rh = sha(root / 'receipt.json')
    assert receipt['status'] == 'complete_functionally_annotated_site_parsimony_exposure_frame'
    bindings = {str(root / 'receipt.json'): rh}
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest
        bindings[str(root / name)] = digest
    for name in ('frame', 'functional', 'projection', 'inputs'):
        path = args[name] / 'receipt.json'
        assert receipt['source_receipts'][name] == sha(path)
        bindings[str(path)] = sha(path)
        source = json.loads(path.read_text())
        for name2, digest in source['artifacts'].items():
            path2 = path.parent / name2
            assert sha(path2) == digest
            bindings[str(path2)] = digest
    for path, digest in receipt['source_readback_sha256'].items():
        assert sha(path) == digest
        bindings[path] = digest
    base = table(args['frame'] / 'site_parsimony_exposure.tsv')
    joined = table(root / 'site_parsimony_exposure_functions.tsv')
    keys = ['marker', 'paired_column_1based']
    assert len(base) == len(joined) == 47529 and not joined.duplicated(keys).any()
    pd.testing.assert_frame_equal(base.sort_values(keys).reset_index(drop=True),
                                  joined[base.columns].sort_values(keys).reset_index(drop=True))
    annotations = table(args['functional'] / 'site_structure_links.tsv')
    observed = annotations[annotations.observed_in_paired_alignment.eq('True')]
    ledger = table(root / 'observed_functional_links.tsv')
    identity = ['marker', 'taxon_id', 'site_id']
    assert not observed.duplicated(identity).any() and not ledger.duplicated(identity).any()
    pd.testing.assert_frame_equal(observed.sort_values(identity).reset_index(drop=True),
                                  ledger[annotations.columns].sort_values(identity).reset_index(drop=True))
    positions = ['marker', 'taxon_id', 'paired_column_1based']
    wanted = set(map(tuple, ledger[positions].to_numpy()))
    found = {}
    with gzip.open(args['projection'] / 'paired_site_accessibility.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[k] for k in positions)
            if key in wanted:
                assert key not in found
                found[key] = row
    assert set(found) == wanted
    for row in ledger.to_dict('records'):
        c = found[tuple(row[k] for k in positions)]
        fields = {'protein_id': 'protein_id', 'model_id': 'model_id',
            'protein_residue_1based': 'protein_residue_1based', 'matrix_column_1based': 'matrix_column_1based',
            'observed_amino_acid': 'amino_acid', 'native_state': '3di_state',
            'coordinate_sasa_angstrom_squared': 'sasa_angstrom_squared',
            'coordinate_ca_plddt': 'ca_plddt', 'coordinate_context': 'context'}
        assert all(row[a] == c[b] for a, b in fields.items())
    groups = {key: part for key, part in observed.groupby(keys)}
    marker_counts = {}
    with_sites = with_candidate = 0
    for row in joined.to_dict('records'):
        key = tuple(row[k] for k in keys)
        part = groups.get(key, observed.iloc[:0])
        taxa = part.taxon_id.nunique()
        candidates = part.loc[part.conserved_candidate.eq('True'), 'taxon_id'].nunique()
        other = part.loc[part.conserved_candidate.eq('False'), 'taxon_id'].nunique()
        expected = dict(functional_correspondence_rows=len(part), functional_correspondence_taxa=taxa,
            conserved_functional_candidate_taxa=candidates, other_functional_correspondence_taxa=other)
        assert all(int(row[k]) == v for k, v in expected.items())
        n = int(row['observed_taxa'])
        assert taxa <= n and n > 0
        assert math.isclose(float(row['functional_correspondence_fraction_observed_taxa']), taxa/n, abs_tol=1e-14)
        assert math.isclose(float(row['conserved_functional_candidate_fraction_observed_taxa']), candidates/n, abs_tol=1e-14)
        assert row['functional_pfam_accessions'] == ';'.join(sorted(set(part.pfam_accession)))
        assert row['functional_annotation_status'] == ('observed_correspondence' if len(part) else 'no_observed_correspondence')
        assert not len(part) or set(part.matrix_column_1based) == {row['matrix_column_1based']}
        with_sites += bool(len(part))
        with_candidate += candidates > 0
        counts = marker_counts.setdefault(row['marker'], Counter())
        counts.update(sites=1, sites_with_correspondence=int(bool(len(part))),
                      sites_with_conserved_candidate=int(candidates > 0), correspondence_rows=len(part))
    summary = table(root / 'marker_summary.tsv')
    assert len(summary) == len(marker_counts) == 125 and summary.marker.is_unique
    for row in summary.to_dict('records'):
        assert all(int(row[k]) == v for k, v in marker_counts[row['marker']].items())
    expected = dict(markers=125, sites=len(base), observed_functional_rows=len(ledger),
        sites_with_correspondence=with_sites, sites_with_conserved_candidate=with_candidate,
        source_site_status_counts=dict(Counter(annotations.site_status)))
    assert all(receipt[k] == value for k, value in expected.items())
    assert sha(opts.plan) == ph
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    result = dict(status='passed_full_recovered_functional_exposure_readback', **expected,
        distinct_annotated_taxon_site_observations=len(wanted), source_receipt_sha256=rh,
        plan_sha256=ph, script_sha256=sha(__file__),
        scope='Every inherited site field, annotation ledger row, coordinate join, per-site count/fraction and marker summary checked. No HMM/ASA recalculation, enrichment, ancestral reconstruction or biological validation.')
    with (root / 'readback.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
