#!/usr/bin/env python3
"""Extract traceable single-copy BUSCO markers; incomplete snapshots are staging only."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
LINEAGE = 'eukaryota_odb12.2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(root, output, allow_incomplete=False):
    with (root / 'metadata/analysis_manifest.tsv').open() as handle:
        taxa = list(csv.DictReader(handle, delimiter='\t'))
    inputs = {r['taxon_id']: r for r in json.loads((root / 'metadata/qc_input_receipts.json').read_text())}
    if len({r['taxon_id'] for r in taxa}) != len(taxa):
        raise ValueError('Duplicate taxon IDs')
    ready, pending = [], []
    for taxon in taxa:
        name = taxon['taxon_id']
        receipt_path = root / f'results/busco/{name}.receipt.json'
        if not receipt_path.exists():
            pending.append(name)
            continue
        receipt = json.loads(receipt_path.read_text())
        if receipt.get('returncode') != 0:
            pending.append(name)
            continue
        source = inputs[name]
        if receipt['input_sha256'] != source['sha256'] or sha(root / source['input_path']) != source['sha256']:
            raise ValueError(f'Stale or changed BUSCO input: {name}')
        ready.append(taxon)
    if pending and not allow_incomplete:
        raise ValueError(f'{len(pending)} taxa await successful QC; use --allow-incomplete only for staging')
    # Publish immutable snapshots. Never silently replace inputs for downstream jobs.
    if output.exists():
        raise FileExistsError(output)
    sequences, mappings, tables = {}, [], []
    expected_markers = None
    for taxon in ready:
        name = taxon['taxon_id']
        run = root / f'results/busco/{name}/run_{LINEAGE}'
        table = run / 'full_table.tsv'
        rows = [line.rstrip('\n').split('\t') for line in table.read_text().splitlines()
                if line and not line.startswith('#')]
        marker_set = {row[0] for row in rows}
        if expected_markers is None:
            expected_markers = marker_set
        if len(marker_set) != 125 or marker_set != expected_markers:
            raise ValueError(f'Inconsistent lineage marker set: {name}')
        with (root / inputs[name]['input_path']).open() as handle:
            source_sequences = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        single = {row[0]: row[2] for row in rows if row[1] == 'Complete'}
        files = {p.stem: p for p in (run / 'busco_sequences/single_copy_busco_sequences').glob('*.faa')}
        if set(single) != set(files):
            raise ValueError(f'Single-copy table/FASTA mismatch: {name}')
        for marker, protein in sorted(single.items()):
            record = SeqIO.read(files[marker], 'fasta')
            seq = str(record.seq)
            if record.id != protein or seq != str(source_sequences[protein].seq):
                raise ValueError(f'Marker differs from source protein: {name}/{marker}')
            sequences.setdefault(marker, {})[name] = seq
            mappings.append({'marker': marker, 'taxon_id': name, 'protein_id': protein,
                             'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                             'source_fasta_sha256': sha(files[marker])})
        tables.append({'taxon_id': name, 'path': str(table.relative_to(root)), 'sha256': sha(table)})
    output.mkdir(parents=True)
    marker_dir = output / 'unaligned'
    marker_dir.mkdir()
    occupancy = []
    roles = {t['taxon_id']: t['study_role'] for t in taxa}
    for marker in sorted(expected_markers or []):
        seqs = sequences.get(marker, {})
        path = marker_dir / f'{marker}.faa'
        path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in sorted(seqs.items())))
        occupancy.append({'marker': marker, 'single_copy_taxa': len(seqs),
                          'ingroup_taxa': sum(roles[n] == 'ingroup' for n in seqs),
                          'outgroup_taxa': sum(roles[n] == 'outgroup' for n in seqs),
                          'fraction_all_planned_taxa': len(seqs) / len(taxa),
                          'sha256': sha(path)})
    for filename, rows, fields in [
        ('protein_mapping.tsv', mappings, ['marker', 'taxon_id', 'protein_id', 'sequence_sha256', 'source_fasta_sha256']),
        ('occupancy.tsv', occupancy, ['marker', 'single_copy_taxa', 'ingroup_taxa', 'outgroup_taxa', 'fraction_all_planned_taxa', 'sha256'])]:
        with (output / filename).open('w') as handle:
            writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    result = {'status': 'incomplete_staging' if pending else 'complete_extraction',
              'planned_taxa': len(taxa), 'completed_taxa': len(ready), 'pending_taxa': pending,
              'lineage': LINEAGE, 'markers': len(occupancy), 'sequences': len(mappings),
              'manifest_sha256': sha(root / 'metadata/analysis_manifest.tsv'),
              'source_tables': tables}
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    return {k: v for k, v in result.items() if k not in ('pending_taxa', 'source_tables')}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    print(json.dumps(prepare(ROOT, args.output, args.allow_incomplete), indent=2))
