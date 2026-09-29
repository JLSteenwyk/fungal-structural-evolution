#!/usr/bin/env python3
"""Bind completed full-cohort geometry, independent readback and RMSD exclusions."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from run_after_verified_dependencies import terminal_state
from screen_duplication_alignment_reuse import sha


def key(row):
    value = row['pair_key'], row['mask'], int(row['order'])
    if value[1] not in ('full', 'plddt70') or value[2] not in (0, 1):
        raise ValueError('Unknown mask/order')
    return value


def launch_plan(launch):
    # Older launch records bind the plan through the recorded argv and digest.
    command = launch['cmdline']
    if command.count('--plan') != 1:
        raise ValueError('Missing or ambiguous launch plan')
    path = command[command.index('--plan')+1]
    if 'plan' in launch and launch['plan'] != path:
        raise ValueError('Conflicting launch plan')
    return path


def indexed(rows):
    result = {}
    for row in rows:
        identity = key(row)
        if identity in result:
            raise ValueError('Duplicate mapping identity')
        result[identity] = (int(row['aligned_length']), row['rmsd_status'])
    return result


def reconcile(geometry_rows, numeric_rows, short_rows, audited_degenerate):
    expected = indexed(numeric_rows)
    shorts = indexed(short_rows)
    audited = indexed(audited_degenerate)
    if any(n not in (1, 2) for n, _ in shorts.values()):
        raise ValueError('Nonshort row in analytic short census')
    counts, rmsd = Counter(), Counter()
    seen, degenerate, observed_shorts = set(), {}, {}
    long_degenerate = []
    for row in geometry_rows:
        identity = key(row)
        item = int(row['aligned_length']), row['rmsd_status']
        if identity in seen or expected.get(identity) != item:
            raise ValueError('Repeated, missing or altered diagnostic mapping')
        seen.add(identity)
        length, classification = item
        if length < 1 or classification not in ('within_printed_rounding', 'outside_printed_rounding'):
            raise ValueError('Unknown RMSD status or empty successful alignment')
        status = row['geometry_status']
        if status not in ('unique_at_numeric_tolerance', 'degenerate_at_numeric_tolerance'):
            raise ValueError('Unknown geometry status')
        counts[identity[1]+':'+status] += 1
        rmsd[classification] += 1
        if status == 'degenerate_at_numeric_tolerance':
            degenerate[identity] = item
            if length >= 3:
                long_degenerate.append(dict(pair_key=identity[0], mask=identity[1],
                    order=identity[2], aligned_length=length, rmsd_status=classification))
        if length < 3:
            if status != 'degenerate_at_numeric_tolerance':
                raise ValueError('Short mapping incorrectly claims unique rotation')
            observed_shorts[identity] = item
    if seen != set(expected) or observed_shorts != shorts or degenerate != audited:
        raise ValueError('Incomplete geometry, short census or audited degeneracy')
    return dict(alignments=len(seen), counts=dict(counts),
        rmsd_classification_counts=dict(rmsd),
        numerically_unique_rotations=len(seen)-len(degenerate),
        degenerate_short_alignments=len(shorts),
        degenerate_long_alignments=len(long_degenerate),
        longer_degenerate_mappings=long_degenerate)


def collect(plan_path):
    plan_path = Path(plan_path)
    plan = json.loads(plan_path.read_text())
    bindings = {str(plan_path): sha(plan_path), **plan['pins']}
    def bind(path, digest):
        path = str(path)
        if path in bindings and bindings[path] != digest:
            raise ValueError('Conflicting source binding: '+path)
        bindings[path] = digest
    def verify():
        for path, digest in bindings.items():
            if sha(path) != digest:
                raise ValueError('Changed source: '+path)
    verify()
    states = {}
    for path in plan['launches']:
        assert path in bindings
        launch = json.loads(Path(path).read_text())
        assert sha(launch_plan(launch)) == launch['plan_sha256']
        state = terminal_state(launch['unit'])
        if state != dict(ActiveState='inactive', Result='success', ExecMainStatus='0'):
            raise ValueError('Unsuccessful geometry dependency: '+launch['unit'])
        states[launch['unit']] = state
    gp = Path(plan['geometry_plan'])
    ap = Path(plan['readback_plan'])
    geometry = json.loads(gp.read_text())
    audit_plan = json.loads(ap.read_text())
    assert audit_plan['source_plan'] == str(gp)
    for spec in (geometry, audit_plan):
        for path, digest in spec['pins'].items():
            bind(path, digest)
    root = Path(geometry['output'])
    receipt_path = root/'receipt.json'
    proof_path = Path(audit_plan['output'])
    diagnostic = Path(geometry['diagnostic'])
    diagnostic_receipt = diagnostic/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    proof = json.loads(proof_path.read_text())
    dr = json.loads(diagnostic_receipt.read_text())
    completion_path = Path(geometry['completion'])
    completion = json.loads(completion_path.read_text())
    short_path = Path(plan['short'])
    short = json.loads(short_path.read_text())
    assert receipt['status'] == 'complete_primary_diagnostic_geometry_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(gp)
    assert proof['status'] == 'passed_full_primary_diagnostic_geometry_readback'
    assert proof['plan_sha256'] == sha(ap)
    assert proof['producer_receipt_sha256'] == sha(receipt_path)
    assert proof['producer_terminal_state'] == dict(ActiveState='inactive', Result='success', ExecMainStatus='0')
    assert dr['status'] == 'complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance'
    assert dr['scientific_eligibility'] is False
    assert receipt['diagnostic_receipt_sha256'] == sha(diagnostic_receipt)
    assert short['status'] == 'passed_all_short_primary_alignment_analytic_rmsd_checks'
    assert short['source_receipt_sha256'] == sha(diagnostic_receipt)
    assert completion['status'] == 'complete_verified_primary_rmsd_diagnostic_and_short_geometry'
    assert completion['scientific_eligibility'] is False
    for path, digest in completion['source_hashes'].items():
        bind(path, digest)
    for path in (receipt_path, proof_path, diagnostic_receipt, completion_path, short_path):
        bind(path, sha(path))
    for folder, record in ((root, receipt), (diagnostic, dr)):
        for name, digest in record['artifacts'].items():
            path = str(folder/name)
            bind(path, digest)
    verify()
    with (root/'alignment_geometry.tsv').open() as gf, (diagnostic/'numeric_readback.tsv').open() as nf:
        result = reconcile(csv.DictReader(gf, delimiter='\t'),
            csv.DictReader(nf, delimiter='\t'), short['rows'], proof['degenerate_alignments'])
    assert result['alignments'] == receipt['alignments'] == proof['alignments_checked'] == dr['numerically_checked_alignments'] == completion['numeric_alignments']
    assert result['counts'] == receipt['counts'] == proof['counts']
    assert result['rmsd_classification_counts'] == dr['rmsd_status_counts']
    assert result['rmsd_classification_counts'].get('outside_printed_rounding', 0) == completion['rmsd_discrepancies']
    assert result['degenerate_short_alignments'] == short['alignments'] == completion['short_alignments']
    result.update(status='complete_verified_primary_diagnostic_geometry', scientific_eligibility=False,
        maximum_scaled_quaternion_curvature_error=proof['maximum_scaled_quaternion_curvature_error'],
        near_zero_quaternion_gaps=proof['near_zero_quaternion_gaps'], terminal_states=states,
        source_hashes=bindings, script_sha256=sha(__file__),
        scope='Every full-cohort geometry identity/length/RMSD classification reconciled with the complete diagnostic table and independent geometry readback. All one/two-residue mappings match the analytic census. Longer numerical degeneracies are retained explicitly rather than excluded by assumption. All original RMSD exclusions remain; no prediction-accuracy, biological effect or inferential acceptance established.')
    verify()
    for unit, state in states.items():
        assert terminal_state(unit) == state
    with Path(plan['output']).open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    result = collect(args.plan)
    print(json.dumps({k:v for k,v in result.items()
        if k not in ('source_hashes', 'longer_degenerate_mappings')}, indent=2))


if __name__ == '__main__':
    main()
