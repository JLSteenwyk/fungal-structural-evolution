#!/usr/bin/env python3
"""Verify full-ledger joins, numeric versions and rehashed corruption rejection."""
import csv
import hashlib
import json
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, record):
    Path(path).write_text(json.dumps(record) + '\n')


def table(path, rows):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def pair(left, right):
    return hashlib.sha256(json.dumps(sorted([left, right]), separators=(',', ':')).encode()).hexdigest()


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        base, catalog, bridge, references = [root / name for name in ['base', 'catalog', 'bridge', 'references']]
        for folder in [base, catalog, bridge, references]:
            folder.mkdir()
        models = []
        for name, version, sequence in [('M', 6, 'AAA'), ('M', 10, 'BBB'), ('Z', 1, 'CCC')]:
            coordinate = root / f'{name}{version}.cif'
            coordinate.write_text(sequence)
            models.append(dict(model_id=name, version=version, sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                               path=str(coordinate), length=3, mean_ca_plddt=90.0))
        for folder, values in [(catalog, models), (base, models[:2])]:
            (folder / 'models.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in values))
        save(catalog / 'receipt.json', dict(artifacts={'models.jsonl': sha(catalog / 'models.jsonl')}))
        table(base / 'model_pairs.tsv', [dict(pair_key=pair(('M', 6), ('M', 10)),
                                            model_a='M', version_a=6, model_b='M', version_b=10)])
        table(base / 'event_model_pair_links.tsv', [dict(guide='profile', gene_a='T0_a', gene_b='T0_b')])
        save(base / 'receipt.json', dict(status='complete_reviewed_duplication_model_pair_queue',
                                        artifacts={p.name: sha(p) for p in base.iterdir()}))
        base_proof = root / 'base-proof.json'
        save(base_proof, dict(status='passed_full_duplication_model_pair_queue_export_readback',
                             source_sha256={str(p): sha(p) for p in base.iterdir()},
                             models_checked=2, distinct_pairs_checked=1))
        database = bridge / 'structure_family_bridge.sqlite'
        with sqlite3.connect(database) as connection:
            connection.execute('CREATE TABLE structures(taxon_id,protein_id,sequence_sha256,model_id,version,model_path)')
            for gene, model in [('T0_a', models[0]), ('T0_b', models[1]), ('T1_r', models[0]), ('T2_r', models[2])]:
                taxon, protein = gene.split('_', 1)
                connection.execute('INSERT INTO structures VALUES(?,?,?,?,?,?)',
                                   (taxon, protein, model['sequence_sha256'], model['model_id'], model['version'], model['path']))
        save(bridge / 'receipt.json', dict(catalog_receipt_sha256=sha(catalog / 'receipt.json')))
        plan_path = root / 'reference-plan.json'
        save(plan_path, dict(bridge=str(database), pins={str(database): sha(database)}))
        source = dict(family='OG', gene_node='n0', gene_a='T0_a', gene_b='T0_b', taxon_id='T0',
                      model_a='M', model_b='M', status='provisional_reference_available',
                      nearest_reference_genes=json.dumps(['T1_r', 'T2_r']), chosen_reference_gene='T1_r',
                      reference_model='M', reference_version=6)
        for guide in ['profile', 'mafft']:
            table(references / (guide + '_sister_references.tsv'), [source])
        guides = [dict(guide=g, candidates=1) for g in ['profile', 'mafft']]
        save(references / 'receipt.json', dict(status='complete_duplication_sister_reference_inventory',
                                              plan_sha256=sha(plan_path), guides=guides,
                                              artifacts={p.name: sha(p) for p in references.iterdir()}))
        reference_proof = root / 'reference-proof.json'
        save(reference_proof, dict(status='passed_full_duplication_sister_reference_readback',
                                  producer_receipt_sha256=sha(references / 'receipt.json'), guides=guides))
        inventory = root / 'inventory'
        subprocess.run([sys.executable, 'scripts/prepare_duplication_reference_comparisons.py',
                        '--inventory', str(references), '--inventory-plan', str(plan_path),
                        '--base-queue', str(base), '--catalog', str(catalog), '--output', str(inventory)],
                       check=True, capture_output=True)
        command = [sys.executable, 'scripts/readback_duplication_reference_comparisons_v2.py',
                   '--inventory', str(inventory), '--references', str(references),
                   '--reference-readback', str(reference_proof), '--base-queue', str(base),
                   '--base-readback', str(base_proof), '--bridge', str(database), '--catalog', str(catalog)]
        proof_path = root / 'proof.json'
        subprocess.run(command + ['--output', str(proof_path)], check=True, capture_output=True)
        proof = json.loads(proof_path.read_text())
        assert proof['event_reference_comparisons'] == 8 and proof['distinct_model_pairs'] == 3
        assert proof['native_gene_model_mappings_checked'] == 4
        assert proof['work_dispositions'] == {'identical_model': 2, 'existing_duplicate_pair': 2, 'additional_pair': 4}
        original = {p.name: p.read_bytes() for p in inventory.iterdir()}

        def reject(number, filename, records):
            table(inventory / filename, records)
            receipt = json.loads((inventory / 'receipt.json').read_text())
            receipt['artifacts'][filename] = sha(inventory / filename)
            save(inventory / 'receipt.json', receipt)
            output = root / f'bad-{number}.json'
            result = subprocess.run(command + ['--output', str(output)], capture_output=True)
            assert result.returncode != 0 and not output.exists()
            for name, data in original.items():
                (inventory / name).write_bytes(data)

        rows = list(csv.DictReader((inventory / 'event_reference_comparisons.tsv').open(), delimiter='\t'))
        # The wrong native gene/model assignment still uses valid catalog models,
        # a valid pair digest and a valid work disposition. It must be rejected.
        rows[0].update(focal_model='Z', focal_version='1', pair_key=pair(('Z', 1), ('M', 6)), work_disposition='additional_pair')
        reject(1, 'event_reference_comparisons.tsv', rows)
        pairs = list(csv.DictReader((inventory / 'model_pairs.tsv').open(), delimiter='\t'))
        pairs[0]['work_disposition'] = 'additional_pair' if pairs[0]['work_disposition'] != 'additional_pair' else 'existing_duplicate_pair'
        reject(2, 'model_pairs.tsv', pairs)
        rows = list(csv.DictReader((inventory / 'event_reference_comparisons.tsv').open(), delimiter='\t'))
        reject(3, 'event_reference_comparisons.tsv', rows[:-1])
        exported = [json.loads(line) for line in (inventory / 'models.jsonl').read_text().splitlines()]
        exported[0]['mean_ca_plddt'] = 10.0
        (inventory / 'models.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in exported))
        receipt = json.loads((inventory / 'receipt.json').read_text())
        receipt['artifacts']['models.jsonl'] = sha(inventory / 'models.jsonl')
        save(inventory / 'receipt.json', receipt)
        output = root / 'bad-metadata.json'
        assert subprocess.run(command + ['--output', str(output)], capture_output=True).returncode != 0
        assert not output.exists()
    print('Complete tied-reference/side grid, numeric versions 6/10, identical/existing/new work, native gene-model joins, and four rehashed corruption rejections passed.')


if __name__ == '__main__':
    main()
