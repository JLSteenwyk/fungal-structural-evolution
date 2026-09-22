#!/usr/bin/env python3
"""Map a completed audited local cohort and independently read every mapped residue."""
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
from audit_busco_gene_copies import ROOT, sha, read_table


def checked(folder):
    receipt = json.loads((folder / 'receipt.json').read_text())
    for name, expected in receipt['artifacts'].items():
        if sha(folder / name) != expected:
            raise ValueError('Changed artifact: ' + name)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_sha = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_sha:
            raise ValueError('Plan changed')
        for name, expected in plan['pins'].items():
            if sha(ROOT / name) != expected:
                raise ValueError('Changed dependency: ' + name)

    verify()
    audit_path, conversion_path = [ROOT / plan[k] for k in ['audit', 'conversion']]
    audit, conversion = checked(audit_path), checked(conversion_path)
    if (audit['status'] != 'complete_artifact_readback' or audit['partial_prediction_snapshot']
            or audit['remaining_eligible'] != 0 or audit['predictions'] != plan['expected_models']
            or conversion['status'] != 'complete_audited_local_model_conversion'
            or conversion['models'] != plan['expected_models']
            or conversion['source_audit_receipt_sha256'] != sha(audit_path / 'receipt.json')):
        raise ValueError('Complete bound audit and conversion required')
    ids = {r['sequence_id'] for r in read_table(audit_path / 'predictions.tsv')}
    inventory = conversion_path / 'inventory.jsonl'
    records = [json.loads(line) for line in inventory.read_text().splitlines()]
    if len(ids) != len(records) or {r['record_id'] for r in records} != ids:
        raise ValueError('Converted model universe differs')
    fields = ['marker', 'taxon_id', 'protein_id', 'sequence_sha256']
    key = lambda r: tuple(r[k] for k in fields)
    global_links = [r for r in read_table(ROOT / plan['global_links']) if 'S' + r['sequence_sha256'] in ids]
    original = read_table(audit_path / 'taxon_links.tsv')
    if Counter(map(key, original)) - Counter(map(key, global_links)) or len(global_links) != plan['expected_global_links']:
        raise ValueError('Global exact-sequence link scope differs')
    control, mapping, readback = [ROOT / plan[k] for k in ['control', 'mapping', 'readback']]
    if any(p.exists() for p in [control, mapping, readback]):
        raise FileExistsError('Use fresh immutable mapping outputs')
    if shutil.disk_usage(ROOT).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient disk headroom')
    control.mkdir(parents=True)
    expected = control / 'expected_global_links.tsv'
    with expected.open('w') as handle:
        writer = csv.DictWriter(handle, list(global_links[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(global_links)
    commands = [
        [sys.executable, 'scripts/map_marker_structures.py', '--inventory', str(inventory),
         '--provider', 'local', '--tool', 'ESMFold v1', '--output', str(mapping)],
        [sys.executable, 'scripts/audit_local_residue_mapping_global_links.py', '--mapping', str(mapping),
         '--prediction-audit', str(audit_path), '--expected-links', str(expected), '--output', str(readback)]]
    (control / 'launch.json').write_text(json.dumps(dict(commands=commands, plan_sha256=plan_sha,
        started_unix=time.time(), expected_global_links_sha256=sha(expected), resources=plan['resources']), indent=2) + '\n')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
    for stage, command in zip(['mapping', 'residue_readback'], commands):
        verify()
        print('Starting', stage, flush=True)
        with (control / (stage + '.log')).open('w') as log:
            subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    mr, rr = checked(mapping), checked(readback)
    if (mr['distinct_models'] != plan['expected_models'] or mr['source_inventory_sha256'] != sha(inventory)
            or rr['status'] != 'passed_full_local_residue_mapping_readback'
            or rr['models'] != plan['expected_models'] or rr['marker_links'] != len(global_links)
            or rr['mapping_receipt_sha256'] != sha(mapping / 'receipt.json')):
        raise ValueError('Mapping or independent readback scope differs')
    result = dict(status='complete_local_cohort_mapping_and_residue_readback',
                  models=rr['models'], marker_links=rr['marker_links'], matrix_residue_links=rr['matrix_residue_links'],
                  plan_sha256=plan_sha, mapping_receipt_sha256=sha(mapping / 'receipt.json'),
                  readback_receipt_sha256=sha(readback / 'receipt.json'),
                  scope='Exact-sequence mapping and every exported retained residue checked. PAE binding, native structural features and evolutionary inference remain separate.')
    (control / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
