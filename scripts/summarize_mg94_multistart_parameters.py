"""Retain parameter ranges across all starts and likelihood-near-best starts."""
import argparse
import csv
import json
import math
import subprocess
import time
from pathlib import Path
import psutil
from audit_local_branch_parameter_profiles import declarations
from readback_whole_proteome_catalog import sha


def parameter_ranges(rows, parameters, tolerance):
    assert len(rows) == len(parameters) and rows
    best = max(r['log_likelihood'] for r in rows)
    near = [i for i, r in enumerate(rows) if best - r['log_likelihood'] <= tolerance]
    names = set(parameters[0])
    assert all(set(p) == names for p in parameters)
    result = []
    for name in sorted(names):
        all_values = [p[name] for p in parameters]
        assert all(math.isfinite(v) and v >= 0 for v in all_values)
        values = [all_values[i] for i in near]
        result.append(dict(parameter=name, all_start_min=min(all_values), all_start_max=max(all_values),
                           near_best_min=min(values), near_best_max=max(values),
                           near_best_range=max(values)-min(values), near_best_starts=len(near)))
    return result


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); digest = sha(args.plan)
    plan = json.loads(args.plan.read_text())
    def verify():
        assert sha(args.plan) == digest
        for path, expected in plan['pins'].items():
            assert sha(path) == expected
    verify(); launch = json.loads(Path(plan['audit_launch']).read_text())
    while True:
        try:
            proc = psutil.Process(launch['pid'])
            live = proc.create_time() == launch['created'] and proc.status() != psutil.STATUS_ZOMBIE
            if live: assert proc.cmdline() == launch['cmdline']
        except psutil.NoSuchProcess:
            live = False
        raw = subprocess.check_output(['systemctl', '--user', 'show', launch['unit'], '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True)
        state = dict(line.split('=', 1) for line in raw.splitlines())
        if state['ActiveState'] == 'failed': raise RuntimeError('Source audit failed')
        if not live and state['ActiveState'] == 'inactive':
            assert state['Result'] == 'success' and state['ExecMainStatus'] == '0'
            break
        print('Waiting for complete codon multistart audit', flush=True); time.sleep(30)
    verify()
    root = Path(plan['producer']); audit_root = Path(plan['audit'])
    receipt = json.loads((root / 'receipt.json').read_text())
    audit = json.loads((audit_root / 'receipt.json').read_text())
    assert audit['status'] == 'passed_full_local_mg94_multistart_artifact_and_numeric_audit'
    assert audit['source_receipt_sha256'] == sha(root / 'receipt.json')
    summary_path = audit_root / 'case_summary.tsv'
    assert sha(summary_path) == audit['artifacts'][summary_path.name]
    with summary_path.open() as h:
        summaries = {r['case_id']: r for r in csv.DictReader(h, delimiter='\t')}
    assert len(summaries) == len(receipt['case_receipts']) == 1632
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False)
    table = out / 'all_parameter_ranges.tsv'
    fields = ['case_id', 'parameter', 'all_start_min', 'all_start_max', 'near_best_min', 'near_best_max', 'near_best_range', 'near_best_starts']
    total = 0; seen = set()
    with table.open('w') as h:
        writer = csv.DictWriter(h, fields, delimiter='\t'); writer.writeheader()
        for proof in receipt['case_receipts']:
            case = proof['case_id']; assert case not in seen; seen.add(case)
            path = root / case / 'receipt.json'; assert sha(path) == proof['receipt_sha256']
            record = json.loads(path.read_text()); rows = record['rows']; assert len(rows) == 8
            proofs = {p['start_label']: p for p in record['proofs']}
            parameters = []
            for row in rows:
                fit = root / case / row['start_label'] / 'fit.bf'
                assert sha(fit) == proofs[row['start_label']]['artifacts']['fit.bf']
                values, fixed = declarations(fit.read_text()); assert not fixed
                parameters.append(values)
            ranges = parameter_ranges(rows, parameters, plan['likelihood_tolerance'])
            assert ranges and ranges[0]['near_best_starts'] == int(summaries[case]['starts_within_1e_5_of_best'])
            assert max(r['log_likelihood'] for r in rows) == float(summaries[case]['best_log_likelihood'])
            for row in ranges:
                writer.writerow(dict(case_id=case, **row)); total += 1
    assert seen == set(summaries)
    # Every serialized range is checked for arithmetic and enclosing bounds.
    count = 0
    with table.open() as h:
        for row in csv.DictReader(h, delimiter='\t'):
            a, b, lo, hi, span = [float(row[k]) for k in fields[2:7]]
            assert a <= lo <= hi <= b and span == hi-lo
            assert 1 <= int(row['near_best_starts']) <= 8
            count += 1
    assert count == total
    verify()
    result = dict(status='complete_full_codon_multistart_parameter_ranges', cases=len(seen), parameter_rows=total,
                  likelihood_tolerance=plan['likelihood_tolerance'], plan_sha256=digest,
                  source_receipt_sha256=sha(root / 'receipt.json'), source_audit_sha256=sha(audit_root / 'receipt.json'),
                  script_sha256=sha(__file__), artifacts={table.name:sha(table)},
                  scope='All free parameters retained across eight starts and within1e-5 likelihood of the best observed start. Descriptive numerical ranges only; not confidence intervals, proof of global optima, identifiability, absence of saturation or selection evidence.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': main()
