"""Unmodified DNA diagnostics under fixed inherited/taxonomic code choices."""
import hashlib
from Bio.Seq import Seq


def digest(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def fixed_roles(target, context):
    old = target['translation_evidence']['translation_table']
    roles = dict(inherited=int(old) if old not in (None, '') else None,
                 snapshot_nuclear=context['nuclear']['code'] or None,
                 snapshot_mitochondrial=context['mitochondrial']['code'] or None)
    assert all(code is None or code in [1, 3, 4, 5, 6, 12, 16, 26] for code in roles.values())
    return roles


def translation(dna, protein, code):
    if len(dna) % 3:
        return dict(status='not_translated_non_triplet_length', code=code,
                    translated_length=None, translated_sha256=None, terminal_stop_removed=None,
                    internal_stop_count=None), None
    raw = str(Seq(dna).translate(table=code))
    stop = raw.endswith('*')
    compared = raw[:-1] if stop else raw
    status = 'translated_no_linked_normalized_protein' if protein is None else (
        'translated_exact_match_to_normalized_protein' if compared == protein else
        'translated_mismatch_to_normalized_protein')
    return dict(status=status, code=code, translated_length=len(compared),
                translated_sha256=digest(compared), terminal_stop_removed=stop,
                internal_stop_count=compared.count('*')), raw


def diagnose(target, dna, protein, context):
    assert set(dna.upper()) <= set('ACGTRYSWKMBDHVN')
    assert digest(dna) == target['original_target_dna_sha256']
    roles = fixed_roles(target, context)
    results, raw = {}, {}
    for code in sorted({c for c in roles.values() if c is not None}):
        results[str(code)], raw[code] = translation(dna, protein, code)
    differences = []
    inherited = roles['inherited']
    if inherited is not None and raw[inherited] is not None:
        for role in ['snapshot_nuclear', 'snapshot_mitochondrial']:
            alternative = roles[role]
            if alternative is None or alternative == inherited:
                continue
            assert len(raw[inherited]) == len(raw[alternative]) == len(dna) // 3
            for offset, (old, new) in enumerate(zip(raw[inherited], raw[alternative])):
                if old != new:
                    differences.append(dict(role=role, inherited_code=inherited, alternative_code=alternative,
                        original_codon_index_1based=offset + 1, original_dna_start_1based=3 * offset + 1,
                        codon=dna[3 * offset:3 * offset + 3], inherited_amino_acid=old,
                        alternative_amino_acid=new,
                        within_normalized_protein_bounds=protein is not None and offset < len(protein),
                        biological_mapping_admitted=False))
    row = dict(ordinal=target['ordinal'], cds_id=target['cds_id'], protein_id=target['protein_id'],
        original_target_dna_sha256=target['original_target_dna_sha256'], original_dna_length=len(dna),
        linked_protein_sha256=digest(protein) if protein is not None else None,
        linked_protein_length=len(protein) if protein is not None else None,
        original_joined_status=target['joined_status'], inherited_translation_evidence=target['translation_evidence'],
        predeclared_code_roles=roles, translation_diagnostics=results,
        changed_codon_counts={role: sum(d['role'] == role for d in differences)
                             for role in ['snapshot_nuclear', 'snapshot_mitochondrial']},
        dna_modified=False, protein_modified=False, best_code_selected=False,
        scientific_eligibility=False, biological_codon_eligibility=False)
    return row, differences
