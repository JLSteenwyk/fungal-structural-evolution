"""Check full analytic-assessment accounting and export remaining review cases.

This is output integrity and independent boundary/flag arithmetic, not a second
evaluation of the likelihood derivative or a proof of global optimality.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
import psutil


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def check_row(row, payload):
    """Return explicit unresolved reasons, retaining every original fit error."""
    if row['original_status'] != payload['status']:
        raise ValueError('Original status mismatch')
    assessment = row['assessment']
    if payload['status'] == 'fit_error_requires_review':
        if assessment != 'original_fit_error_unassessed' or row['error_type'] != payload['error_type']:
            raise ValueError('Fit error lost or changed')
        return ['original_fit_error']
    if payload['status'] not in ('candidate_passed_numerical_optimization_checks',
                                 'candidate_requires_optimization_review'):
        raise ValueError('Unknown production status')
    if assessment == 'analytic_diagnostic_unresolved':
        if not row.get('error_type') or not isinstance(row.get('error'), str):
            raise ValueError('Missing unresolved diagnostic')
        return ['analytic_diagnostic_unresolved']
    if assessment != 'analytic_score_evaluated':
        raise ValueError('Unknown assessment')
    theta = payload['log1p_ratios']
    gradient = row['analytic_gradient']
    reported = row['analytic_projected_gradient']
    if not all(len(v) == 3 for v in (theta, gradient, reported)):
        raise ValueError('Wrong gradient dimensions')
    if not all(math.isfinite(x) for v in (theta, gradient, reported) for x in v):
        raise ValueError('Nonfinite gradient')
    upper = math.log1p(payload['maximum_ratio'])
    expected = []
    for position, derivative in zip(theta, gradient):
        if not 0 <= position <= upper:
            raise ValueError('Parameter outside bounds')
        if position <= 1e-7:
            derivative = min(derivative, 0.)
        elif position >= upper - 1e-7:
            derivative = max(derivative, 0.)
        expected.append(derivative)
    if reported != expected:
        raise ValueError('Boundary projection differs')
    passed = max(map(abs, expected)) <= 1e-3
    flags = payload['checks']
    other = bool(flags['best_optimizer_success'] and not flags['upper_bound_contact']
                 and flags['all_full_face_starts_agree'])
    checks = {'analytic_threshold_pass': passed,
              'other_optimizer_checks_pass': other,
              'original_gradient_pass': flags['projected_gradient_pass'],
              'all_numerical_checks_with_analytic_gradient': passed and other}
    if any(type(row[k]) is not bool or row[k] != value for k, value in checks.items()):
        raise ValueError('Assessment flag mismatch')
    reasons = []
    if not passed:
        reasons.append('analytic_gradient_threshold')
    if not flags['best_optimizer_success']:
        reasons.append('optimizer_termination')
    if flags['upper_bound_contact']:
        reasons.append('upper_bound_contact')
    if not flags['all_full_face_starts_agree']:
        reasons.append('full_face_start_disagreement')
    return reasons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_hash:
            raise ValueError('Plan changed')
        for path, expected in config['pins'].items():
            if sha(path) != expected:
                raise ValueError('Pinned file changed: ' + path)

    verify()
    launch = json.loads(Path(config['launch']).read_text())
    if launch['plan_sha256'] != sha(config['assessment_plan']):
        raise ValueError('Launch/plan mismatch')
    while True:
        try:
            process = psutil.Process(launch['pid'])
            if process.create_time() != launch['created'] or process.status() == psutil.STATUS_ZOMBIE:
                break
            if process.cmdline() != launch['cmdline']:
                raise ValueError('Dependency command changed')
        except psutil.NoSuchProcess:
            break
        time.sleep(30)
    raw = subprocess.check_output(['systemctl', '--user', 'show', launch['unit'],
                                   '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True)
    state = dict(line.split('=', 1) for line in raw.splitlines())
    if state != {'ActiveState': 'inactive', 'Result': 'success', 'ExecMainStatus': '0'}:
        raise ValueError('Assessment not successfully terminal')
    verify()
    assessment_plan = json.loads(Path(config['assessment_plan']).read_text())
    for path, expected in assessment_plan['pins'].items():
        if sha(path) != expected:
            raise ValueError('Assessment input changed: ' + path)
    production_plan = json.loads(Path(assessment_plan['production_plan']).read_text())
    root = Path(production_plan['output'])
    source = Path(assessment_plan['output'])
    receipt = json.loads((source / 'receipt.json').read_text())
    if (receipt['status'] != 'complete_full_analytic_stationarity_assessment_pending_readback'
            or receipt['plan_sha256'] != sha(config['assessment_plan'])
            or receipt['production_receipt_sha256'] != sha(root / 'receipt.json')
            or receipt['production_audit_sha256'] != sha(assessment_plan['audit'])):
        raise ValueError('Assessment receipt binding failed')
    audit = json.loads(Path(assessment_plan['audit']).read_text())
    if (audit['status'] != 'passed_full_working_model_output_integrity_audit'
            or audit['source_receipt_sha256'] != sha(root / 'receipt.json')):
        raise ValueError('Production audit binding failed')
    production = json.loads((root / 'receipt.json').read_text())
    manifest_path = root / 'fit_manifest.jsonl'
    if sha(manifest_path) != production['artifacts']['fit_manifest.jsonl']:
        raise ValueError('Changed production manifest')
    manifest = {}
    for line in manifest_path.open():
        entry = json.loads(line)
        key = entry['fit_input_id'], entry['tree']
        if key in manifest:
            raise ValueError('Duplicate production key')
        manifest[key] = entry
    recipes = [json.loads(line)['fit_input_id'] for line in
               (Path(production_plan['inventory']) / 'unique_fit_recipes.jsonl').open()]
    trees = {p.stem for p in Path(production_plan['factors']).glob('*.npz')}
    expected = {(identifier, tree) for identifier in recipes for tree in trees}
    if set(manifest) != expected or len(expected) != 144040:
        raise ValueError('Incomplete production scope')
    data = source / 'analytic_stationarity.jsonl'
    if sha(data) != receipt['artifacts'][data.name]:
        raise ValueError('Changed analytic table')
    out = Path(config['output'])
    out.mkdir(parents=True, exist_ok=False)
    seen = set()
    counts, transitions, reasons_count = Counter(), Counter(), Counter()
    review_count = 0
    with (out / 'remaining_review_cases.jsonl').open('x') as queue:
        for line in data.open():
            row = json.loads(line)
            key = row['fit_input_id'], row['tree']
            if key not in manifest or key in seen:
                raise ValueError('Extra or duplicate assessment')
            seen.add(key)
            entry = manifest[key]
            if row['source_fit_sha256'] != entry['sha256'] or sha(entry['path']) != entry['sha256']:
                raise ValueError('Source fit changed')
            saved = json.loads(Path(entry['path']).read_text())
            if saved['fit_input_id'] != key[0] or saved['tree'] != key[1]:
                raise ValueError('Fit identity mismatch')
            reasons = check_row(row, saved['payload'])
            counts[row['assessment']] += 1
            if row['assessment'] == 'analytic_score_evaluated':
                transitions[str(row['original_gradient_pass']) + '->' + str(row['analytic_threshold_pass'])] += 1
            if reasons:
                review_count += 1
                reasons_count.update(reasons)
                queue.write(json.dumps(dict(fit_input_id=key[0], tree=key[1],
                    source_fit=entry['path'], source_fit_sha256=entry['sha256'],
                    reasons=reasons), allow_nan=False) + '\n')
            if len(seen) % 10000 == 0:
                print('Checked analytic dispositions', len(seen), '/ 144040', flush=True)
    if (seen != expected or receipt['dispositions'] != len(seen)
            or dict(counts) != receipt['assessment_counts']
            or dict(transitions) != receipt['gradient_threshold_transitions']):
        raise ValueError('Assessment totals or full coverage mismatch')
    if sha(data) != receipt['artifacts'][data.name]:
        raise ValueError('Analytic table changed during readback')
    verify()
    result = dict(status='passed_full_analytic_stationarity_output_readback',
        dispositions=len(seen), assessment_counts=dict(counts),
        gradient_threshold_transitions=dict(transitions), remaining_review_cases=review_count,
        review_reason_counts=dict(reasons_count), source_receipt_sha256=sha(source / 'receipt.json'),
        plan_sha256=plan_hash, script_sha256=sha(__file__),
        artifacts={'remaining_review_cases.jsonl': sha(out / 'remaining_review_cases.jsonl')},
        scope='Every assessment and immutable production source checked; scalar boundary projection and all flags independently reconstructed. Remaining review cases retained without editing originals. Does not recompute derivatives, prove global optima, validate uncertainty or establish biological effects.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
