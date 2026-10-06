#!/usr/bin/env python3
"""Exhaustive independent codon semantics and unmodified DNA controls."""
import argparse
from datetime import datetime, timezone
from itertools import product
import json
from pathlib import Path
import Bio
import Bio.Data.CodonTable
import Bio.Data.IUPACData
from Bio.Seq import Seq

from independent_codon_translation_v1 import DNA, lookup, translate
from unmodified_cds_translation_v1 import diagnose, digest, translation
from reference_measurement_union_sources import bind


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receipt', type=Path, required=True)
    args = p.parse_args()
    assert not args.receipt.exists()
    context_path = Path('metadata/full_taxonomic_genetic_code_context_20261006_v1.json')
    codebook = json.loads(context_path.read_text())['codebook']
    checks = 0
    for code in [1, 3, 4, 5, 6, 12, 16, 26]:
        table = lookup(codebook, code)
        for symbols in product(DNA, repeat=3):
            dna = ''.join(symbols)
            assert translate(dna, table) == str(Seq(dna).translate(table=code)), (code, dna)
            checks += 1
    assert checks == 27000
    controls = []
    for dna, protein, code, expected in [
        ('ATGCTGTAA', 'ML', 1, 'translated_exact_match_to_normalized_protein'),
        ('ATGCTGTAA', 'ML', 12, 'translated_mismatch_to_normalized_protein'),
        ('GTGTAA', 'M', 1, 'translated_mismatch_to_normalized_protein'),
        ('ATGTAATAA', 'M', 1, 'translated_mismatch_to_normalized_protein'),
        ('ATGT', 'M', 1, 'not_translated_non_triplet_length'),
        ('ATGTAA', None, 1, 'translated_no_linked_normalized_protein')]:
        result, raw = translation(dna, protein, code)
        assert result['status'] == expected
        if dna == 'ATGTAATAA':
            assert result['internal_stop_count'] == 1 and result['translated_length'] == 2
        controls.append(dna + ':' + expected)
    dna = 'ATGCTGTAA'
    target = dict(ordinal=1, cds_id='literal', protein_id='protein', original_target_dna_sha256=digest(dna),
                  joined_status='original_exact', translation_evidence=dict(translation_table='1',
                  code_provenance='table_1_assumption_no_explicit_code', inherited_status='exact_translation'))
    context = dict(nuclear=dict(code=12), mitochondrial=dict(code=4))
    row, sites = diagnose(target, dna, 'ML', context)
    assert row['predeclared_code_roles'] == dict(inherited=1, snapshot_nuclear=12, snapshot_mitochondrial=4)
    assert row['changed_codon_counts']['snapshot_nuclear'] == 1
    assert sites == [dict(role='snapshot_nuclear', inherited_code=1, alternative_code=12,
                         original_codon_index_1based=2, original_dna_start_1based=4, codon='CTG',
                         inherited_amino_acid='L', alternative_amino_acid='S',
                         within_normalized_protein_bounds=True, biological_mapping_admitted=False)]
    assert not row['dna_modified'] and not row['protein_modified'] and not row['best_code_selected']
    controls.append('fixed_code_comparison_preserves_original_exact_without_selecting_best_code')
    dna2 = 'ATGTAA'
    row2, sites2 = diagnose(dict(target, original_target_dna_sha256=digest(dna2)), dna2, 'M',
                           dict(nuclear=dict(code=6), mitochondrial=dict(code=0)))
    assert sites2[0]['original_codon_index_1based'] == 2 and not sites2[0]['within_normalized_protein_bounds']
    assert row2['predeclared_code_roles']['snapshot_mitochondrial'] is None
    controls.append('terminal_codon_difference_not_mapped_to_nonexistent_protein_residue')
    try:
        diagnose(dict(target, original_target_dna_sha256='wrong'), dna, 'ML', context)
    except AssertionError:
        controls.append('reject_changed_original_dna')
    else:
        raise AssertionError('source DNA digest accepted')
    pins = {}
    for q in [Path(__file__), Path('scripts/unmodified_cds_translation_v1.py'),
              Path('scripts/independent_codon_translation_v1.py'), context_path,
              Path(Bio.__file__), Path(__import__('Bio.Seq', fromlist=['__file__']).__file__),
              Path(Bio.Data.CodonTable.__file__), Path(Bio.Data.IUPACData.__file__)]:
        bind(pins, q)
    result = dict(status='passed_exhaustive_fixed_code_translation_and_unmodified_source_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), exhaustive_IUPAC_codon_code_checks=checks,
        fixed_codes=[1, 3, 4, 5, 6, 12, 16, 26], controls=controls, control_count=len(controls),
        biopython_version=Bio.__version__, source_hashes=pins, scientific_eligibility=False,
        biological_codon_eligibility=False, gpu=False, new_predictions=0,
        scope='All3375IUPAC codons in eight fixed source/taxonomic codes compared between installed '
              'Biopython and independent pinned-codebook expansion. Literal source/hash/terminal/partial '
              'and no-code-selection controls. Shared original codebook/library dependencies explicit; '
              'no real taxon pilot or biological code/compartment/selection admission.')
    with args.receipt.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':
    main()
