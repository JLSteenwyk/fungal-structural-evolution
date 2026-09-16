#!/usr/bin/env python3
"""Verify every full-domain query link against source representative FASTA files."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def records(path):
    """Independent minimal FASTA parser; reject sequence before a header."""
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                if name is None or any(c.isspace() for c in line):
                    raise ValueError('Invalid FASTA sequence line')
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def query_ids(path):
    ids, residues = set(), 0
    for name, sequence in records(path):
        h = hashlib.sha256(sequence.encode()).digest()
        require(sequence and name == 'S' + h.hex() and h not in ids,
                'Empty, duplicate, or misidentified query')
        ids.add(h)
        residues += len(sequence)
    return ids, residues


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    require(not args.output.exists(), 'Use a fresh output directory')
    started = time.time()
    root = args.root
    full = root / 'data/domains/full-inputs-v1'
    marker = root / 'data/domains/marker-inputs-v1'
    source = root / 'results/gene_representatives/full-v2/receipt.json'
    manifest = root / 'metadata/analysis_manifest.tsv'
    receipt = json.loads((full / 'receipt.json').read_text())
    sr = json.loads(source.read_text())
    mr = json.loads((marker / 'receipt.json').read_text())
    pins = {str(p.relative_to(root)): digest(p) for p in
            [full / 'receipt.json', marker / 'receipt.json', source, manifest]}
    for path, expected in [(source, receipt['source_receipt_sha256']),
                           (manifest, receipt['manifest_sha256']),
                           (marker / 'receipt.json', receipt['marker_input_receipt_sha256']),
                           (marker / 'sequences.faa', mr['artifacts']['sequences.faa'])]:
        require(digest(path) == expected, 'Changed source: ' + str(path))
    for name, expected in receipt['artifacts'].items():
        require(digest(full / name) == expected, 'Changed full-domain input')
    with manifest.open() as handle:
        taxa = [x['taxon_id'] for x in csv.DictReader(handle, delimiter='\t')]
    require(sr['status'] == 'complete' and len(taxa) == len(set(taxa))
            and set(taxa) == {x['taxon_id'] for x in sr['taxa']}, 'Taxon universe differs')
    additional, residues = query_ids(full / 'additional_sequences.faa')
    markers, _ = query_ids(marker / 'sequences.faa')
    require(not additional.intersection(markers), 'Search partitions overlap')
    require(len(additional) == receipt['additional_unique_sequences']
            and residues == receipt['additional_unique_residues'], 'Query totals differ')
    print('Query sequences verified', len(additional), len(markers), flush=True)
    seen_additional, used_markers = set(), set()
    summaries = []
    expected_taxa = {x['taxon_id']: x for x in receipt['taxon_inputs']}
    require(set(expected_taxa) == set(taxa), 'Receipt taxon universe differs')
    with (full / 'protein_links.tsv').open() as handle:
        links = iter(csv.DictReader(handle, delimiter='\t'))
        for item in sorted(sr['taxa'], key=lambda x: x['taxon_id']):
            taxon = item['taxon_id']
            path = root / item['path']
            require(digest(path) == item['sha256'] == expected_taxa[taxon]['source_sha256'],
                    'Changed representative FASTA')
            proteins, reused = set(), 0
            for protein, sequence in records(path):
                require(sequence and protein not in proteins, 'Empty/duplicate protein')
                proteins.add(protein)
                h = hashlib.sha256(sequence.encode()).digest()
                partition = 'marker-inputs-v1' if h in markers else 'additional_full_proteome'
                if h in markers:
                    reused += 1
                    used_markers.add(h)
                else:
                    require(h in additional, 'Source sequence missing from queries')
                    seen_additional.add(h)
                row = next(links, None)
                require(row == dict(taxon_id=taxon, protein_id=protein,
                                    sequence_id='S' + h.hex(), query_source=partition),
                        'Protein link missing, changed, duplicated, or reordered')
            expected = expected_taxa[taxon]
            require(len(proteins) == item['selected_proteins'] == expected['proteins']
                    and reused == expected['marker_query_links']
                    and len(proteins) - reused == expected['additional_query_links'],
                    'Taxon count mismatch')
            summaries.append(dict(taxon_id=taxon, proteins_checked=len(proteins),
                                  marker_query_links=reused,
                                  additional_query_links=len(proteins) - reused))
            if len(summaries) % 50 == 0:
                print('Taxa checked', len(summaries), flush=True)
        require(next(links, None) is None, 'Extra protein links')
    require(seen_additional == additional, 'Orphan additional queries')
    counts = dict(taxa=len(summaries), representative_proteins=sum(x['proteins_checked'] for x in summaries),
                  additional_unique_sequences=len(additional), reused_unique_marker_sequences=len(used_markers),
                  full_proteome_unique_sequences=len(additional) + len(used_markers),
                  marker_query_protein_links=sum(x['marker_query_links'] for x in summaries))
    require(all(receipt[k] == v for k, v in counts.items()), 'Aggregate count mismatch')
    require(all(digest(root / p) == h for p, h in pins.items()), 'Metadata changed during readback')
    args.output.mkdir(parents=True)
    table = args.output / 'taxon_link_readback.tsv'
    with table.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(summaries)
    result = dict(status='passed_full_domain_protein_link_readback', **counts,
                  additional_query_residues_checked=residues,
                  marker_queries_outside_representative_proteomes=len(markers - used_markers),
                  elapsed_seconds=time.time() - started, source_receipt_pins=pins,
                  script_sha256=digest(Path(__file__)), artifacts={table.name: digest(table)},
                  scope='All protein links and query sequence hashes independently reconstructed from representative FASTA files. Exact source and partition coverage checked. Does not validate gene-representative selection, taxonomy, homology, domain assignments, or biological absence.')
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
