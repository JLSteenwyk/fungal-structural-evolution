#!/usr/bin/env python3
"""Verify published ORF coordinates by exact translation against deposited genomes."""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
TAXA = [('F2109901', 'Amoeboradix_gromovi_K1', 'Amoeboradix_gromovi_K1'),
        ('F2020955', 'Sanchytrium_tribonematis_SC1', 'Sanchytrium_tribonematis_SC')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(record, genomes, aliases):
    result = {'protein_id': record.id, 'status': 'unparsed_header', 'contig': '',
              'start': '', 'end': '', 'strand': '', 'contig_name_match': '', 'orf_id': ''}
    match = re.search(r' (\S+):(\d+)-(\d+)\(([+-])\)$', record.description)
    if not match:
        return result, None
    contig, a, b, strand = match.groups()
    low, high = sorted((int(a), int(b)))
    result.update(start=low, end=high, strand=strand)
    if contig in genomes:
        target = contig
        result['contig_name_match'] = 'exact'
    else:
        candidates = aliases.get(contig.split('_cov_')[0], [])
        if len(candidates) != 1:
            result['status'] = 'missing_or_ambiguous_contig'
            return result, None
        target = candidates[0]
        result['contig_name_match'] = 'node_and_length_prefix_only'
    result['contig'] = target
    if low < 1 or high > len(genomes[target]):
        result['status'] = 'coordinates_out_of_bounds'
        return result, None
    sequence = genomes[target].seq[low - 1:high]
    if strand == '-':
        sequence = sequence.reverse_complement()
    if len(sequence) % 3:
        result['status'] = 'not_codon_multiple'
        return result, None
    if str(sequence.translate(table=1)).rstrip('*') != str(record.seq).rstrip('*'):
        result['status'] = 'translation_mismatch'
        return result, None
    result.update(status='exact_translation', orf_id=f'ORF:{target}:{low}-{high}:{strand}')
    return result, str(sequence)


def main():
    sources = {r['path']: r for r in json.loads((ROOT / 'metadata/external_genome_receipts.json').read_text())}
    folder = ROOT / 'results/sanchytrid_coordinates'
    folder.mkdir(parents=True, exist_ok=True)
    summaries = []
    for taxon, protein_name, genome_name in TAXA:
        protein = ROOT / f'data/external/16411629/{protein_name}_proteome.fasta'
        genome = ROOT / f'data/external/14816085/{genome_name}_genome_assembly_contigs.fasta'
        for path in [protein, genome]:
            if sha(path) != sources[str(path.relative_to(ROOT))]['sha256']:
                raise ValueError(f'Changed source: {path}')
        with genome.open() as handle:
            genomes = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        aliases = defaultdict(list)
        for contig in genomes:
            aliases[contig.split('_cov_')[0]].append(contig)
        rows = []
        cds_path = folder / f'{taxon}.verified_cds.fna'
        with protein.open() as handle, cds_path.open('w') as cds:
            for record in SeqIO.parse(handle, 'fasta'):
                row, sequence = verify(record, genomes, aliases)
                rows.append(row)
                if sequence is not None:
                    cds.write(f'>{record.id}\n{sequence}\n')
        path = folder / f'{taxon}.tsv'
        with path.open('w') as out:
            writer = csv.DictWriter(out, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        summary = {'taxon_id': taxon, 'proteins': len(rows), 'status_counts': dict(Counter(r['status'] for r in rows)),
                   'contig_aliases_verified': sum(r['status'] == 'exact_translation' and r['contig_name_match'] != 'exact' for r in rows),
                   'protein_sha256': sha(protein), 'genome_sha256': sha(genome),
                   'mapping_path': str(path.relative_to(ROOT)), 'mapping_sha256': sha(path),
                   'cds_path': str(cds_path.relative_to(ROOT)), 'cds_sha256': sha(cds_path)}
        summaries.append(summary)
        print(taxon, summary['status_counts'], flush=True)
    (ROOT / 'metadata/sanchytrid_coordinate_audit.json').write_text(json.dumps({
        'taxa': summaries, 'translation_table': 1,
        'note': 'Verified genomic ORFs are not resolved gene/isoform models. Overlap, annotation completeness and assembly discrepancies still require assessment.'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
