"""Check rigid-fit identifiability for every strictly audited primary alignment."""
import argparse
from collections import Counter
import csv
from functools import lru_cache
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import psutil
from assess_reference_alignment_geometry import geometry
from readback_reference_alignment_geometry_v2 import verify_row
from duplication_alignment_numeric_readback import load_pdb, check_alignment
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for path, expected in plan['pins'].items():
            if sha(path) != expected:
                raise ValueError('Pinned input changed: ' + path)

    verify()
    launch = json.loads(Path(plan['audit_launch']).read_text())
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != launch['cmdline']:
                raise ValueError('Audit process identity changed')
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    raw = subprocess.check_output(['systemctl', '--user', 'show', launch['service'],
        '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True)
    state = dict(line.split('=', 1) for line in raw.splitlines())
    if state != dict(ActiveState='inactive', Result='success', ExecMainStatus='0'):
        raise ValueError('Strict primary audit did not finish successfully; no bypass permitted')
    verify()
    audit_plan = json.loads(Path(launch['plan']).read_text())
    if sha(launch['plan']) != launch['plan_sha256'] or audit_plan['mode'] != 'primary':
        raise ValueError('Wrong audit plan')
    source_plan = json.loads(Path(audit_plan['source_plan']).read_text())
    root = Path(source_plan['output'])
    audited = Path(audit_plan['output'])
    audit_path = audited / 'receipt.json'
    audit = json.loads(audit_path.read_text())
    producer_path = root / 'receipt.json'
    producer = json.loads(producer_path.read_text())
    if (audit['status'] != 'passed_full_duplication_alignment_mapping_rmsd_identity_readback'
            or audit['mode'] != 'primary' or audit['plan_sha256'] != sha(launch['plan'])
            or audit['producer_receipt_sha256'] != sha(producer_path)
            or producer['plan_sha256'] != sha(audit_plan['source_plan'])):
        raise ValueError('Unbound strict audit')
    table = audited / 'numeric_readback.tsv'
    manifest_path = root / 'checkpoint_manifest.tsv'
    if (sha(table) != audit['artifacts'][table.name]
            or sha(manifest_path) != producer['artifacts'][manifest_path.name]):
        raise ValueError('Changed source table')
    manifest = {}
    for item in csv.DictReader(manifest_path.open(), delimiter='\t'):
        if item['path'] in manifest:
            raise ValueError('Duplicate checkpoint')
        manifest[item['path']] = item
    expected = {path for path, item in manifest.items() if item['status'] == 'aligned'}
    if len(manifest) != audit['directed_dispositions'] or len(expected) != audit['numerically_checked_alignments']:
        raise ValueError('Incomplete audit scope')
    folder = Path(source_plan['inputs'])
    input_receipt = json.loads((folder / 'receipt.json').read_text())
    input_path = folder / 'inputs.jsonl'
    if (sha(input_path) != input_receipt['artifacts'][input_path.name]
            or input_receipt['plan_sha256'] != sha(source_plan['input_plan'])):
        raise ValueError('Changed materialized inputs')
    inputs = {}
    for line in input_path.open():
        row = json.loads(line)
        key = row['model_id'], row['version'], row['mask']
        if key in inputs:
            raise ValueError('Duplicate input model')
        inputs[key] = row

    @lru_cache(maxsize=128)
    def coordinates(key):
        return load_pdb(inputs[key])

    bindings = {str(p): sha(p) for p in (audit_path, producer_path, table,
                manifest_path, input_path, folder / 'receipt.json')}
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    output_table = out / 'alignment_geometry.tsv'
    seen = set()
    counts = Counter()
    maximum_error = 0.
    near_zero = 0
    fields = ['pair_key', 'mask', 'order'] + list(geometry(np.zeros((3, 3)), np.zeros((3, 3))))
    fields += ['quaternion_scaled_curvature_error', 'near_zero_quaternion_gap']
    with output_table.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in csv.DictReader(table.open(), delimiter='\t'):
            pair, mask, order = row['pair_key'], row['mask'], int(row['order'])
            rel = f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json'
            if rel not in expected or rel in seen:
                raise ValueError('Extra or duplicate numeric row')
            seen.add(rel)
            if sha(root / rel) != manifest[rel]['sha256']:
                raise ValueError('Changed native alignment')
            record = json.loads((root / rel).read_text())
            if (record['pair_key'], record['mask'], record['order'], record['status']) != (pair, mask, order, 'aligned'):
                raise ValueError('Wrong native alignment identity')
            left, right = [coordinates((r['model_id'], r['version'], mask)) for r in record['inputs']]
            numeric = check_alignment(record, left, right)
            if any(float(row[name]) != float(value) for name, value in numeric.items()):
                raise ValueError('Numeric audit mismatch')
            i = j = 0
            ix, iy = [], []
            for a, b in zip(record['metrics']['alignment_left'], record['metrics']['alignment_right']):
                if a != '-' and b != '-':
                    ix.append(i)
                    iy.append(j)
                i += a != '-'
                j += b != '-'
            x, y = left[1][ix], right[1][iy]
            metrics = geometry(x, y)
            error, close = verify_row(metrics, x, y)
            counts[mask + ':' + metrics['geometry_status']] += 1
            maximum_error = max(maximum_error, error)
            near_zero += int(close)
            writer.writerow(dict(pair_key=pair, mask=mask, order=order, **metrics,
                quaternion_scaled_curvature_error=error, near_zero_quaternion_gap=close))
            if len(seen) % 10000 == 0:
                print('Primary geometry checked', len(seen), '/', len(expected), flush=True)
    if seen != expected:
        raise ValueError('Missing primary geometry rows')
    for path, expected_hash in bindings.items():
        if sha(path) != expected_hash:
            raise ValueError('Source changed during geometry assessment')
    verify()
    receipt = dict(status='complete_primary_alignment_geometry_with_quaternion_checks',
        plan_sha256=plan_hash, alignments=len(seen), directed_dispositions=audit['directed_dispositions'],
        skipped_dispositions=audit['directed_dispositions'] - len(seen), counts=dict(counts),
        maximum_scaled_quaternion_curvature_error=maximum_error, near_zero_quaternion_gaps=near_zero,
        source_sha256=bindings, artifacts={output_table.name: sha(output_table)},
        scope='All strictly audited successful primary mappings rechecked from hashed PDBs. SVD geometry checked with alternate SVD and quaternion eigensystem for every alignment. Numerical identifiability only, not prediction uncertainty, structural asymmetry, model adequacy or biological inference. Skipped native dispositions remain upstream. Serialized table readback remains pending.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
