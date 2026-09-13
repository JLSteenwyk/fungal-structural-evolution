#!/usr/bin/env python3
"""Prepare full within-genus/code diagnostic alignments with explicit caveats."""
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['codons', 'coverage', 'hybrids', 'marker_review', 'output']:
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable diagnostic input directory')
    cr, gr = checked_receipt(a.codons), checked_receipt(a.coverage)
    if gr['source_codon_receipt_sha256'] != sha(a.codons / 'receipt.json'):
        raise ValueError('Codon/coverage snapshots differ')
    hybrids = {r['taxon_id'] for r in read_table(a.hybrids)}
    review = {}
    for r in read_table(a.marker_review):
        marker, flag = r['marker'], r['curated_orthology_caveat']
        if marker in review and review[marker] != flag:
            raise ValueError('Conflicting marker caveats')
        review[marker] = flag
    source_columns = json.loads((a.codons / 'protein_source_columns.json').read_text())
    source_audit = {(r['marker'], r['taxon_id']): r for r in read_table(a.codons / 'sequence_audit.tsv') if r['status'] == 'translation_verified_codon_alignment'}
    groups = [r for r in read_table(a.coverage / 'group_codon_coverage.tsv') if r['policy'] == 'exclude_recorded_annotation_gene_label_flags']
    a.output.mkdir(parents=True); summaries = []; exclusions = []; memberships = []; seen = set()
    cache = {}; total_observed = 0
    for group in groups:
        marker, genus, code = group['marker'], group['genus_label'], group['translation_table']
        case = genus + '__' + marker + '__code' + code
        if case in seen or not all(c.isalnum() or c in '_-' for c in case):
            raise ValueError('Duplicate or unsafe case identity')
        seen.add(case)
        original = json.loads(group['retained_taxa_json'])
        if len(original) != len(set(original)):
            raise ValueError('Duplicate source taxon')
        taxa = sorted(set(original) - hybrids)
        for taxon in original:
            memberships.append({'case_id': case, 'taxon_id': taxon, 'retained': taxon in taxa,
                                'reason': 'retained_prior_strict_policy' if taxon in taxa else 'curated_hybrid_requires_homeolog_analysis'})
        if marker not in cache:
            records = list(SeqIO.parse(a.codons / (marker + '.fna'), 'fasta'))
            cache[marker] = {r.id: str(r.seq) for r in records}
            if len(cache[marker]) != len(records):
                raise ValueError('Duplicate source codon sequence')
        sequences = {t: cache[marker][t] for t in taxa}
        width = len(source_columns[marker]); kept = []; encoded = {}; counts = {}
        if taxa:
            table = CodonTable.unambiguous_dna_by_id[int(code)]
            for taxon, sequence in sequences.items():
                sr = source_audit[marker, taxon]
                if sr['translation_table'] != code or len(sequence) != width*3:
                    raise ValueError('Code or source dimensions differ')
                codons = [sequence[i:i+3] for i in range(0,len(sequence),3)]
                if any(set(c) <= set('ACGT') and c not in table.forward_table for c in codons):
                    raise ValueError('Unexpected stop codon')
                encoded[taxon] = codons
            for i in range(width):
                called = sum(set(encoded[t][i]) <= set('ACGT') for t in taxa)
                if called >= math.ceil(.8*len(taxa)):
                    kept.append(i)
            counts = {t: sum(set(encoded[t][i]) <= set('ACGT') for i in kept) for t in taxa}
        ready = len(taxa) >= 4 and len(kept) >= 100 and all(counts.values())
        row = {'case_id': case, 'marker': marker, 'genus_label': genus, 'translation_table': code,
               'source_strict_taxa': len(original), 'excluded_curated_hybrids': len(set(original)&hybrids),
               'retained_taxa': len(taxa), 'source_codon_columns': width, 'retained_codon_columns': len(kept),
               'minimum_observed_codons': min(counts.values()) if counts else 0,
               'maximum_observed_codons': max(counts.values()) if counts else 0,
               'marker_copy_caveat': review.get(marker,'none_recorded'),
               'status': 'ready_for_tree_and_divergence_diagnostics' if ready else 'insufficient_coverage_after_taxon_policy',
               'selection_status': 'not_validated_for_selection_testing'}
        summaries.append(row)
        if not ready:
            exclusions.append(row); continue
        folder = a.output / case; folder.mkdir()
        with (folder/'codons.fna').open('w') as dna, (folder/'amino_acids.faa').open('w') as protein:
            for taxon in taxa:
                codons = [encoded[taxon][i] if set(encoded[taxon][i]) <= set('ACGT') else '???' for i in kept]
                aa = ''.join(table.forward_table[c] if c != '???' else '?' for c in codons)
                dna.write('>'+taxon+'\n'+''.join(codons)+'\n'); protein.write('>'+taxon+'\n'+aa+'\n')
                total_observed += counts[taxon]
        write_table(folder/'columns.tsv',[{'diagnostic_codon_column_1based':j+1,'source_codon_column_1based':i+1,'source_mafft_protein_column_1based':source_columns[marker][i]} for j,i in enumerate(kept)])
    write_table(a.output/'case_summary.tsv',summaries)
    write_table(a.output/'membership.tsv',memberships)
    if exclusions:
        write_table(a.output/'excluded_cases.tsv',exclusions)
    result = {'status':'complete_genus_codon_diagnostic_inputs','cases_audited':len(summaries),
              'ready_cases':sum(r['status']=='ready_for_tree_and_divergence_diagnostics' for r in summaries),
              'observed_taxon_codons_in_ready_cases':total_observed,
              'status_counts':dict(Counter(r['status'] for r in summaries)),
              'source_receipts':{'codons':sha(a.codons/'receipt.json'),'coverage':sha(a.coverage/'receipt.json')},
              'hybrid_table_sha256':sha(a.hybrids),'marker_review_sha256':sha(a.marker_review),
              'script_sha256':sha(Path(__file__)),
              'interpretation':'All previously strict genus/code groups retained in the ledger; curated hybrids excluded from ordinary gene-tree diagnostics. Complete called codons retained at >=80% taxon occupancy, invariant sites retained, unknown codons represented as ???. Ready requires >=4 taxa, >=100 columns and nonzero observed codons in each taxon. Source MAFFT correspondence preserved, not new alignment validation. Genus names do not establish monophyly, copy identity, independent transitions, absence of saturation or selection suitability. TFIIB copy caveat retained. These are diagnostic inputs, not selection results.',
              'artifacts':{str(p.relative_to(a.output)):sha(p) for p in sorted(a.output.rglob('*')) if p.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2),flush=True)


if __name__ == '__main__':
    main()
