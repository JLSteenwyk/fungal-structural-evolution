#!/usr/bin/env python3
"""Close full common-core validation and count order-robust coverage by model multiplicity."""
import csv
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
import psutil
from screen_duplication_domain_alignment_coverage import sha


def main():
    pp = Path('metadata/whole_protein_common_fits_plan_20260927.json')
    plan = json.loads(pp.read_text())
    dependency = json.loads(Path('metadata/whole_protein_common_fits_readback_launch_20260927.json').read_text())
    while True:
        try:
            proc = psutil.Process(dependency['pid'])
            if proc.create_time() != dependency['created'] or proc.status() == psutil.STATUS_ZOMBIE:
                break
            assert proc.cmdline() == dependency['cmdline']
        except psutil.NoSuchProcess:
            break
        print('Waiting for full common-core audit', dependency['pid'], flush=True)
        time.sleep(30)
    states = {}
    sources = [pp]
    for suffix in ['launch', 'readback_launch']:
        path = Path(f'metadata/whole_protein_common_fits_{suffix}_20260927.json')
        launch = json.loads(path.read_text()); sources.append(path)
        assert launch['plan_sha256'] == sha(pp)
        try:
            proc = psutil.Process(launch['pid'])
            assert proc.create_time() != launch['created'], 'Original process still present'
        except psutil.NoSuchProcess:
            pass
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState',
             '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), state
        states[launch['unit']] = state
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    root = Path(plan['output']); rp = root/'receipt.json'
    ap = Path('metadata/whole_protein_common_fits_readback_20260927.json')
    r = json.loads(rp.read_text()); a = json.loads(ap.read_text()); sources += [rp, ap]
    assert r['status'] == 'complete_common_residue_whole_protein_fits_pending_readback'
    assert a['status'] == 'passed_full_common_residue_whole_protein_fit_readback'
    assert r['plan_sha256'] == a['plan_sha256'] == sha(pp)
    assert a['producer_receipt_sha256'] == sha(rp)
    assert r['fit_rows'] == a['fit_rows_checked'] == 563808
    assert r['mapping_dispositions'] == a['mapping_dispositions'] == 281904
    for key in ['counts', 'screen_pass_counts']:
        assert r[key] == a[key]
    for name, digest in r['artifacts'].items():
        assert sha(root/name) == digest
    table = root/'common_residue_fits.tsv'
    screens = [s['id'] for s in plan['screens']]
    keys = ['triad_id', 'distinct_models', 'mask', 'mapping_definition']
    frame = pd.read_csv(table, sep='\t', usecols=keys+[s+'_pass' for s in screens])
    grouped = frame.groupby(keys, sort=True)
    assert (grouped.size() == 8).all() and len(grouped) == 70476
    rows = []
    for screen in screens:
        totals = grouped[screen+'_pass'].sum().reset_index(name='passing_orders')
        for key, group in totals.groupby(keys[1:], sort=True):
            rows.append(dict(zip(keys[1:], key), screen=screen,
                             total_triads=len(group),
                             all_eight_orders_pass=int((group.passing_orders == 8).sum()),
                             some_orders_pass=int(group.passing_orders.between(1, 7).sum()),
                             no_orders_pass=int((group.passing_orders == 0).sum())))
    # Independently aggregate the entire serialized table without pandas.
    scalar = defaultdict(Counter); counts = Counter()
    with table.open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = tuple(row[k] for k in keys)
            counts[key] += 1
            for screen in screens:
                value = int(row[screen+'_pass']); assert value in [0, 1]
                scalar[key][screen] += value
    assert len(counts) == 70476 and set(counts.values()) == {8}
    expected = defaultdict(Counter)
    for key, values in scalar.items():
        for screen, count in values.items():
            bucket = 'all_eight_orders_pass' if count == 8 else 'no_orders_pass' if count == 0 else 'some_orders_pass'
            expected[key[1:] + (screen,)][bucket] += 1
            expected[key[1:] + (screen,)]['total_triads'] += 1
    for row in rows:
        key = tuple(str(row[k]) for k in keys[1:]) + (row['screen'],)
        for field in ['total_triads', 'all_eight_orders_pass', 'some_orders_pass', 'no_orders_pass']:
            assert row[field] == expected[key][field], (key, field)
    assert len(rows) == 72
    out = Path('results/structural_comparisons/whole-protein-common-fit-summary-20260927-v1')
    out.mkdir(exist_ok=False)
    target = out/'order_robust_coverage.tsv'
    with target.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    with target.open() as f:
        actual = list(csv.DictReader(f, delimiter='\t'))
    assert actual == [{k: str(v) for k, v in row.items()} for row in rows]
    assert sha(table) == r['artifacts']['common_residue_fits.tsv']
    receipt = dict(status='complete_verified_whole_protein_common_fits_and_coverage_summary',
                   terminal_states=states, source_hashes={str(p): sha(p) for p in sources},
                   script_sha256=sha(__file__), fit_rows=563808, triad_mask_definition_groups=70476,
                   counts=a['counts'], screen_pass_counts=a['screen_pass_counts'],
                   maximum_absolute_rmsd_or_contrast_difference=a['maximum_absolute_rmsd_or_contrast_difference'],
                   summary_rows=len(rows), artifacts={str(target): sha(target)},
                   scope='Full independent quaternion readback and two complete coverage aggregations. Counts are oriented model triads, not independent duplication events. Separate 1/2/3-model strata; all orders retained. Coverage robustness is not evidence of biological asymmetry or an ancestral reference.')
    with Path('metadata/whole_protein_common_fits_completed_20260927.json').open('x') as f:
        json.dump(receipt, f, indent=2); f.write('\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
