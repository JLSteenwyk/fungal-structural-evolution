#!/usr/bin/env python3
"""Map exact protein accessions through GFF Parent links; flag ambiguous loci."""
import csv
import fcntl
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def attributes(text):
    result = {}
    for field in text.split(';'):
        if not field or field == '.':
            continue
        key, value = field.split('=', 1)
        # Split delimiters before decoding escaped literal commas.
        result[key] = [unquote(item) for item in value.split(',')]
    return result


def parse_gff(path, transcript_ids=False):
    parents = defaultdict(set)
    genes, proteins = set(), defaultdict(set)
    partial = set()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            attrs = attributes(fields[8])
            if transcript_ids and fields[2] not in ('gene', 'pseudogene', 'mRNA', 'transcript'):
                continue
            ids = attrs.get('ID', [])
            if len(ids) != 1:
                if fields[2] in ('gene', 'pseudogene', 'CDS'):
                    raise ValueError('Relevant feature lacks a unique ID')
                continue
            ident = ids[0]
            parents[ident].update(attrs.get('Parent', []))
            if fields[2] in ('gene', 'pseudogene'):
                genes.add(ident)
            if transcript_ids and fields[2] in ('mRNA', 'transcript'):
                proteins[ident].add(ident)
            if fields[2] == 'CDS':
                for protein in attrs.get('protein_id', []):
                    proteins[protein].add(ident)
                    if attrs.get('partial') == ['true']:
                        partial.add(protein)

    def ancestors(ident, trail):
        if ident in trail:
            raise ValueError('Cycle in GFF Parent links')
        if ident in genes:
            return {ident}
        found = set()
        for parent in parents.get(ident, set()):
            found.update(ancestors(parent, trail | {ident}))
        return found

    mapping = {protein: set().union(*(ancestors(ident, set()) for ident in ids))
               for protein, ids in proteins.items()}
    return mapping, partial


def parse_creolimax_gtf(path):
    mapping = defaultdict(set)
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'CDS':
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', fields[8]))
            if not attrs.get('gene_id') or not attrs.get('transcript_id'):
                raise ValueError('CDS lacks gene/transcript ID')
            mapping[attrs['transcript_id']].add(attrs['gene_id'])
    return dict(mapping), set()


def main():
    folder = ROOT / 'results/gene_mapping'
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder / '.mapping.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    with (ROOT / 'metadata/analysis_manifest.tsv').open() as handle:
        taxa = {row['taxon_id'] for row in csv.DictReader(handle, delimiter='\t')}
    inputs = {r['taxon_id']: r for r in json.loads((ROOT / 'metadata/qc_input_receipts.json').read_text())}
    annotations = {}
    for line in (ROOT / 'data/raw/annotation_downloads.jsonl').read_text().splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue  # A concurrent writer may be appending the final line.
        if record['status'] == 'validated' and record['taxon_id'] in taxa:
            annotations[record['taxon_id']] = record
    external = {r['path']: r for r in json.loads((ROOT / 'metadata/external_genome_receipts.json').read_text())}
    for article, basename in [('5426470', 'Clim_long.annot.gff'), ('5426494', 'Nk52_long.annot.gff'),
                              ('5426506', 'Pgem_long.annot.gff'), ('5426458', 'Awhi_long.annot.gff'),
                              ('1403592', 'Creolimax_fragrantissima.gtf.gz')]:
        name = 'OFS' + article
        path = f'data/external/{article}/{basename}'
        if name in taxa:
            annotations[name] = dict(external[path], taxon_id=name,
                                     mapping_mode='creolimax_gtf' if article == '1403592' else 'transcript_gff')
    coordinate_audit = ROOT / 'metadata/sanchytrid_coordinate_audit.json'
    if coordinate_audit.exists():
        for record in json.loads(coordinate_audit.read_text())['taxa']:
            if record['taxon_id'] in taxa:
                annotations[record['taxon_id']] = dict(record, path=record['mapping_path'],
                    sha256=record['mapping_sha256'], mapping_mode='verified_orf_coordinates')
    summaries = []
    for name, annotation in sorted(annotations.items()):
        gff = ROOT / annotation['path']
        source = ROOT / inputs[name]['input_path']
        for path, expected in [(gff, annotation['sha256']), (source, inputs[name]['sha256'])]:
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f'Changed source: {path}')
        mode = annotation.get('mapping_mode', 'ncbi_protein_gff')
        if mode == 'verified_orf_coordinates':
            with gff.open() as handle:
                mapping = {r['protein_id']: {r['orf_id']} for r in csv.DictReader(handle, delimiter='\t')
                           if r['status'] == 'exact_translation'}
            partial = set()
        else:
            mapping, partial = parse_creolimax_gtf(gff) if mode == 'creolimax_gtf' else parse_gff(gff, mode == 'transcript_gff')
        rows, counts = [], Counter()
        with source.open() as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                loci = mapping.get(record.id, set())
                status = 'unique_gene' if len(loci) == 1 else ('multiple_genes' if loci else 'unmapped')
                if mode == 'verified_orf_coordinates' and loci:
                    status = 'provisional_orf'
                counts[status] += 1
                rows.append({'taxon_id': name, 'protein_id': record.id,
                             'gene_ids_json': json.dumps(sorted(loci)), 'status': status,
                             'partial_cds': 'unknown' if mode in ('verified_orf_coordinates', 'transcript_gff', 'creolimax_gtf') else record.id in partial,
                             'protein_length': len(record.seq)})
        target = folder / f'{name}.tsv'
        with target.with_suffix('.partial').open('w') as out:
            writer = csv.DictWriter(out, list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        target.with_suffix('.partial').replace(target)
        gene_counts = Counter(json.loads(row['gene_ids_json'])[0] for row in rows if row['status'] == 'unique_gene')
        summaries.append({'taxon_id': name, 'proteins': len(rows), **dict(counts),
                          'uniquely_mapped_genes': len(gene_counts),
                          'genes_with_multiple_proteins': sum(n > 1 for n in gene_counts.values()),
                          'mapping_mode': mode,
                          'annotation_sha256': annotation['sha256'], 'proteome_sha256': inputs[name]['sha256'],
                          'mapping_path': str(target.relative_to(ROOT)),
                          'mapping_sha256': hashlib.sha256(target.read_bytes()).hexdigest()})
    (ROOT / 'metadata/gene_mapping_snapshot.json').write_text(json.dumps({
        'planned_taxa': len(taxa), 'mapped_taxa': len(summaries), 'taxa': summaries,
        'note': 'No isoform selection performed; ambiguous/unmapped proteins require resolution. Gene IDs are local to this annotation version.'}, indent=2) + '\n')
    print('Mapped', len(summaries), 'taxa; multi-protein genes:', sum(r['genes_with_multiple_proteins'] for r in summaries))


if __name__ == '__main__':
    main()
