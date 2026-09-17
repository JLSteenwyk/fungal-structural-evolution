#!/usr/bin/env python3
"""Check every staged sequence against original source FASTAs with Bio.SeqIO."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        assert sha(Path(path)) == digest, path
    root = Path(plan['output'])
    receipt = json.loads((root / 'receipt.json').read_text())
    assert receipt['plan_sha256'] == sha(a.plan)
    table = root / 'family_sequences.tsv'
    assert sha(table) == receipt['artifacts'][table.name]
    with (Path(plan['census']) / 'unique_tree_tasks.tsv').open() as handle:
        tasks = {r['membership_sha256']: r for r in csv.DictReader(handle, delimiter='\t')
                 if r['source_type'] == 'discovery'}
    with table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert {r['membership_sha256'] for r in rows} == set(tasks)
    assert len(rows) == len(tasks) == receipt['families']
    assert {p.relative_to(root).as_posix() for p in (root / 'families').iterdir()} == {r['relative_path'] for r in rows}
    observed = {}
    records = residues = total_bytes = 0
    for row in rows:
        key = row['membership_sha256']
        path = root / row['relative_path']
        assert not path.is_symlink() and sha(path) == row['sha256']
        sequences = list(SeqIO.parse(path, 'fasta'))
        names = [r.id for r in sequences]
        assert len(set(names)) == len(names) == int(row['proteins']) == int(tasks[key]['proteins'])
        assert hashlib.sha256(('\n'.join(sorted(names)) + '\n').encode()).hexdigest() == key
        assert len({n.split('_', 1)[0] for n in names}) == int(row['taxa'])
        lengths = [len(r.seq) for r in sequences]
        assert sum(lengths) == int(row['residues'])
        assert min(lengths) == int(row['minimum_length']) and max(lengths) == int(row['maximum_length'])
        assert path.stat().st_size == int(row['bytes'])
        for sequence in sequences:
            value = hashlib.sha256(str(sequence.seq).encode()).hexdigest()
            if sequence.id in observed:
                assert observed[sequence.id] == value
            observed[sequence.id] = value
        records += len(names)
        residues += sum(lengths)
        total_bytes += path.stat().st_size
    with (Path(plan['source_root']) / 'copied_files.tsv').open() as handle:
        sources = [r for r in csv.DictReader(handle, delimiter='\t')
                   if Path(r['relative_path']).name.startswith('Species') and r['relative_path'].endswith('.fa')]
    matched = set()
    scanned = 0
    for row in sources:
        path = Path(row['source'])
        assert sha(path) == row['sha256'], path
        for record in SeqIO.parse(path, 'fasta'):
            scanned += 1
            if record.id in observed:
                assert record.id not in matched
                assert hashlib.sha256(str(record.seq).encode()).hexdigest() == observed[record.id]
                matched.add(record.id)
    assert matched == set(observed)
    assert len(matched) == receipt['distinct_proteins']
    assert scanned == receipt['source_proteins_scanned'] == 5815847
    assert records == receipt['family_sequence_records']
    assert residues == receipt['residues'] and total_bytes == receipt['fasta_bytes']
    result = dict(status='passed_complete_new_family_sequence_readback',
                  families=len(rows), distinct_proteins=len(matched), family_sequence_records=records,
                  residues=residues, original_source_files=len(sources),
                  staging_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(Path(__file__)),
                  scope='All staged native IDs and full sequences independently checked against original source FASTAs; exact task membership, dimensions and file hashes checked. No alignment or tree validation.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
