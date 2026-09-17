#!/usr/bin/env python3
"""Stage complete native-ID sequences for every distinct newly discovered tree task."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from census_expanded_family_tree_workload import clusters, sha


def fasta(path):
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    for path, digest in plan['pins'].items():
        if sha(Path(path)) != digest:
            raise ValueError('Changed pinned source: ' + path)
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free < plan['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient free disk')
    census = Path(plan['census'])
    cr = json.loads((census / 'receipt.json').read_text())
    ca = json.loads((census / 'readback.json').read_text())
    if ca['status'] != 'passed_native_expanded_tree_workload_readback' or ca['census_receipt_sha256'] != sha(census / 'receipt.json'):
        raise ValueError('Unverified census')
    if sha(census / 'unique_tree_tasks.tsv') != cr['artifacts']['unique_tree_tasks.tsv']:
        raise ValueError('Changed task list')
    with (census / 'unique_tree_tasks.tsv').open() as handle:
        tasks = [r for r in csv.DictReader(handle, delimiter='\t') if r['source_type'] == 'discovery']
    if len(tasks) != plan['expected_tasks']:
        raise ValueError('Unexpected number of tasks')
    merged = Path(plan['merged'])
    mr = json.loads((merged / 'receipt.json').read_text())
    wanted = {(t['representative_guide'], t['representative_family']): t for t in tasks}
    genes_by_task = {}
    needed = set()
    for guide in mr['guides']:
        path = merged / guide['guide'] / 'clusters_id_pairs.txt'
        if sha(path) != guide['artifacts'][str(path.relative_to(merged))]:
            raise ValueError('Changed partition')
        for family, genes in clusters(path):
            task = wanted.get((guide['guide'], family))
            if task is None:
                continue
            key = task['membership_sha256']
            if hashlib.sha256(('\n'.join(sorted(genes)) + '\n').encode()).hexdigest() != key:
                raise ValueError('Task membership differs')
            if len(genes) != int(task['proteins']):
                raise ValueError('Task size differs')
            genes_by_task[key] = sorted(genes)
            needed.update(genes)
    if len(genes_by_task) != len(tasks):
        raise ValueError('Missing task memberships')
    source_root = Path(plan['source_root'])
    with (source_root / 'copied_files.tsv').open() as handle:
        sources = [r for r in csv.DictReader(handle, delimiter='\t')
                   if Path(r['relative_path']).name.startswith('Species') and r['relative_path'].endswith('.fa')]
    if len(sources) != 526:
        raise ValueError('Incomplete source proteomes')
    sequences = {}
    scanned = 0
    for row in sources:
        path = source_root / row['relative_path']
        if sha(path) != row['sha256']:
            raise ValueError('Changed source proteome')
        for name, seq in fasta(path):
            scanned += 1
            if name in needed:
                if name in sequences or not seq:
                    raise ValueError('Duplicate or empty required sequence')
                sequences[name] = seq
    if scanned != 5815847 or set(sequences) != needed:
        raise ValueError('Incomplete source sequence coverage')
    out.mkdir(parents=True)
    (out / 'families').mkdir()
    table = out / 'family_sequences.tsv'
    total_records = total_residues = total_bytes = 0
    with table.open('w') as handle:
        fields = ['membership_sha256', 'relative_path', 'proteins', 'taxa', 'residues',
                  'minimum_length', 'maximum_length', 'sha256', 'bytes']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for task in tasks:
            key = task['membership_sha256']
            path = out / 'families' / (key + '.faa')
            genes = genes_by_task[key]
            lengths = [len(sequences[g]) for g in genes]
            with path.open('w') as output:
                for g in genes:
                    output.write('>' + g + '\n' + sequences[g] + '\n')
            size = path.stat().st_size
            writer.writerow(dict(membership_sha256=key, relative_path=str(path.relative_to(out)),
                                 proteins=len(genes), taxa=task['taxa'], residues=sum(lengths),
                                 minimum_length=min(lengths), maximum_length=max(lengths),
                                 sha256=sha(path), bytes=size))
            total_records += len(genes)
            total_residues += sum(lengths)
            total_bytes += size
            if total_bytes > plan['output_allowance_gib'] * 2**30:
                raise ValueError('Output allowance exceeded')
    result = dict(status='complete_new_family_sequence_staging_requires_readback',
                  families=len(tasks), source_proteins_scanned=scanned,
                  distinct_proteins=len(needed), family_sequence_records=total_records,
                  residues=total_residues, fasta_bytes=total_bytes,
                  plan_sha256=sha(a.plan), script_sha256=sha(Path(__file__)),
                  artifacts={table.name: sha(table)},
                  scope='All distinct new tree tasks; native identifiers and complete source sequences preserved. No sequence deduplication, trimming, alignment or tree inference.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
