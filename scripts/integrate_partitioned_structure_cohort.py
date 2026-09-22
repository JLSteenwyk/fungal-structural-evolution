#!/usr/bin/env python3
"""Join validated disjoint model inventories, then map and independently audit residues."""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import psutil
from Bio import SeqIO
from audit_busco_gene_copies import ROOT, sha, read_table


def read(path):
    return json.loads(path.read_text())


def checked(folder):
    receipt = read(folder / 'receipt.json')
    for name, digest in receipt['artifacts'].items():
        if sha(folder / name) != digest:
            raise ValueError('Changed artifact: ' + str(folder / name))
    return receipt


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan, plan_sha = read(a.plan), sha(a.plan)

    def verify():
        if sha(a.plan) != plan_sha:
            raise ValueError('Plan changed')
        for name, digest in plan['pins'].items():
            if sha(ROOT / name) != digest:
                raise ValueError('Changed dependency: ' + name)

    verify()
    control, union_audit, inventory_dir, mapping, readback = [ROOT / plan[k] for k in
        ['control', 'union_audit', 'inventory', 'mapping', 'readback']]
    if any(p.exists() for p in [control, union_audit, inventory_dir, mapping, readback]):
        raise FileExistsError('Use fresh immutable outputs')
    control.mkdir(parents=True)
    state = dict(status='waiting_for_partition_validation', plan_sha256=plan_sha)

    def save():
        state['updated_unix'] = time.time()
        temporary = control / 'state.tmp'
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(control / 'state.json')

    save()
    try:
        while True:
            try:
                proc = psutil.Process(plan['predecessor_pid'])
                live = proc.create_time() == plan['predecessor_create_time'] and proc.status() != psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:
                live = False
            if not live:
                break
            time.sleep(20)
        verify()
        handoff_path = ROOT / plan['handoff'] / 'receipt.json'
        handoff = read(handoff_path)
        if (handoff['status'] != 'complete_partitioned_prediction_snapshot_handoff'
                or handoff['config_sha256'] != sha(ROOT / plan['handoff_config'])
                or handoff['predictions'] != plan['expected_models']):
            raise ValueError('Partition handoff did not complete cleanly')
        if shutil.disk_usage(ROOT).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
            raise ValueError('Insufficient disk headroom')
        state['status'] = 'combining_validated_partitions'
        save()
        all_predictions, all_links, all_models, source_records = [], [], [], []
        seen = set()
        for cohort in handoff['cohorts']:
            source_audit, source_conversion = ROOT / cohort['audit'], ROOT / cohort['conversion']
            ar, cr = checked(source_audit), checked(source_conversion)
            if (sha(source_audit / 'receipt.json') != cohort['audit_receipt_sha256']
                    or sha(source_conversion / 'receipt.json') != cohort['conversion_receipt_sha256']
                    or cr['source_audit_receipt_sha256'] != cohort['audit_receipt_sha256']
                    or cr['status'] != 'complete_audited_local_model_conversion'
                    or ar['predictions'] != cohort['predictions'] or cr['models'] != cohort['predictions']):
                raise ValueError('Partition provenance or dimensions differ')
            predictions = read_table(source_audit / 'predictions.tsv')
            ids = {r['sequence_id'] for r in predictions}
            models = [json.loads(line) for line in (source_conversion / 'inventory.jsonl').read_text().splitlines()]
            if (len(ids) != len(predictions) or ids & seen or len(ids) != cohort['predictions']
                    or len(models) != len(ids) or {r['record_id'] for r in models} != ids):
                raise ValueError('Overlapping or incomplete partition identity set')
            for record in models:
                if record['status'] != 'verified' or len(record['models']) != 1:
                    raise ValueError('Unexpected converted model record')
                model = record['models'][0]
                if (model['prediction_config_sha256'] != ar['config_sha256']
                        or model['sequence_sha256'] != record['record_id'][1:]
                        or model['provider'] != 'local' or model['tool'] != 'ESMFold v1'):
                    raise ValueError('Prediction configuration or source identity differs')
            seen.update(ids)
            all_predictions.extend(predictions)
            all_links.extend(read_table(source_audit / 'taxon_links.tsv'))
            all_models.extend(models)
            source_records.append(dict(cohort, audit_status=ar['status'],
                                       source_prediction_config_sha256=ar['config_sha256']))
        expected_records = list(SeqIO.parse(ROOT / plan['original_fasta'], 'fasta'))
        expected_ids = {r.id for r in expected_records}
        if seen != expected_ids or len(expected_records) != len(seen) or len(seen) != plan['expected_models']:
            raise ValueError('Combined models do not exactly cover the original queue')
        fields = ['marker', 'taxon_id', 'protein_id', 'sequence_sha256']
        key = lambda row: tuple(row[k] for k in fields)
        expected_original = [r for r in read_table(ROOT / plan['original_links']) if r['sequence_id'] in seen]
        global_links = [r for r in read_table(ROOT / plan['global_links']) if 'S' + r['sequence_sha256'] in seen]
        if (Counter(map(key, all_links)) != Counter(map(key, expected_original))
                or Counter(map(key, all_links)) - Counter(map(key, global_links))
                or len(global_links) != plan['expected_global_links']):
            raise ValueError('Original or global marker-link multiplicity differs')
        union_audit.mkdir(parents=True)
        write_table(union_audit / 'predictions.tsv', sorted(all_predictions, key=lambda r: r['sequence_id']))
        write_table(union_audit / 'taxon_links.tsv', all_links)
        (union_audit / 'partition_sources.json').write_text(json.dumps(source_records, indent=2) + '\n')
        union_receipt = dict(status='complete_validated_partitioned_artifact_union', predictions=len(seen),
                             taxon_marker_links=len(all_links), handoff_receipt_sha256=sha(handoff_path),
                             artifacts={p.name: sha(p) for p in union_audit.iterdir()},
                             scope='Complete disjoint union of previously independently audited partitions; original partial-source status preserved in partition_sources.json. No single shared configuration hash is invented.')
        (union_audit / 'receipt.json').write_text(json.dumps(union_receipt, indent=2) + '\n')
        inventory_dir.mkdir(parents=True)
        inventory = inventory_dir / 'inventory.jsonl'
        inventory.write_text(''.join(json.dumps(r) + '\n' for r in sorted(all_models, key=lambda r: r['record_id'])))
        (inventory_dir / 'receipt.json').write_text(json.dumps(dict(status='complete_verified_partitioned_model_inventory',
            models=len(seen), union_audit_receipt_sha256=sha(union_audit / 'receipt.json'),
            artifacts={'inventory.jsonl': sha(inventory)}), indent=2) + '\n')
        expected_path = control / 'expected_global_links.tsv'
        write_table(expected_path, global_links)
        commands = [
            [sys.executable, 'scripts/map_marker_structures.py', '--inventory', str(inventory),
             '--provider', 'local', '--tool', 'ESMFold v1', '--output', str(mapping)],
            [sys.executable, 'scripts/audit_local_residue_mapping_global_links.py', '--mapping', str(mapping),
             '--prediction-audit', str(union_audit), '--expected-links', str(expected_path), '--output', str(readback)]]
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
        for stage, command in zip(['mapping', 'residue_readback'], commands):
            verify()
            state.update(status='running_stage', stage=stage, command=command)
            save()
            with (control / (stage + '.log')).open('w') as log:
                subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        mr, rr = checked(mapping), checked(readback)
        if (mr['distinct_models'] != len(seen) or rr['status'] != 'passed_full_local_residue_mapping_readback'
                or rr['models'] != len(seen) or rr['marker_links'] != len(global_links)
                or rr['mapping_receipt_sha256'] != sha(mapping / 'receipt.json')):
            raise ValueError('Mapping/readback scope differs')
        verify()
        result = dict(status='complete_local_cohort_mapping_and_residue_readback',
                      models=len(seen), marker_links=len(global_links), matrix_residue_links=rr['matrix_residue_links'],
                      plan_sha256=plan_sha, mapping_receipt_sha256=sha(mapping / 'receipt.json'),
                      readback_receipt_sha256=sha(readback / 'receipt.json'),
                      union_audit_receipt_sha256=sha(union_audit / 'receipt.json'),
                      scope='All partition-specific model provenance retained; complete original queue and global links mapped and independently read back. Native features and confidence qualification remain pending.')
        (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        state['status'] = result['status']
        save()
    except Exception as exc:
        state.update(status='failed_requires_review', error=repr(exc))
        save()
        raise


if __name__ == '__main__':
    main()
