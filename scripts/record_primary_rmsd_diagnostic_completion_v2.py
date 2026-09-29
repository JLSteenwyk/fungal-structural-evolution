#!/usr/bin/env python3
"""Bind full discrepancy accounting and short-mapping checks from a pinned plan."""
import argparse
import json
from pathlib import Path
from run_after_verified_dependencies import terminal_state
from screen_duplication_domain_alignment_coverage import sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    bindings = {str(args.plan): sha(args.plan), **plan['pins']}
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    states = {}
    for path in plan['launches']:
        assert path in bindings
        launch = json.loads(Path(path).read_text())
        if 'plan' in launch:
            assert sha(launch['plan']) == launch['plan_sha256']
        state = terminal_state(launch['unit'])
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), (launch['unit'], state)
        assert launch['unit'] not in states
        states[launch['unit']] = state
    assert len(states) == 3
    root = Path(plan['diagnostic'])
    rp = root/'receipt.json'
    cp = Path(plan['accounting'])
    sp = Path(plan['short'])
    receipt, check, short = [json.loads(p.read_text()) for p in [rp, cp, sp]]
    assert receipt['status'] == 'complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance'
    assert receipt['scientific_eligibility'] is False
    assert check['status'] == 'verified_completed_primary_rmsd_diagnostic_not_scientific_acceptance'
    assert short['status'] == 'passed_all_short_primary_alignment_analytic_rmsd_checks'
    assert check['receipt_sha256'] == short['source_receipt_sha256'] == sha(rp)
    assert receipt['plan_sha256'] == check['plan_sha256'] == sha(plan['diagnostic_plan'])
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    assert check['numeric_rows'] == receipt['numerically_checked_alignments']
    assert check['numeric_rows'] + check['explicit_unavailable'] == receipt['directed_dispositions']
    lengths = check['aligned_length_counts']
    assert short['alignments'] == sum(lengths.get(str(n), 0) for n in [1, 2])
    assert short['one_residue'] == lengths.get('1', 0)
    assert short['two_residue'] == lengths.get('2', 0)
    key = lambda r: (r['pair_key'], r['mask'], r['order'])
    bad = {key(r) for r in receipt['rmsd_discrepancies']}
    assert len(bad) == len(receipt['rmsd_discrepancies'])
    assert len({key(r) for r in short['rows']}) == len(short['rows']) == short['alignments']
    bad_short = {key(r) for r in short['rows'] if r['rmsd_status'] == 'outside_printed_rounding'}
    assert bad_short <= bad and len(bad_short) == short['outside_printed_rounding']
    verify()
    result = dict(status='complete_verified_primary_rmsd_diagnostic_and_short_geometry',
        scientific_eligibility=False, terminal_states=states,
        source_hashes={str(p): sha(p) for p in [rp, cp, sp, args.plan, Path(plan['diagnostic_plan'])]},
        script_sha256=sha(__file__), dispositions=receipt['directed_dispositions'],
        numeric_alignments=check['numeric_rows'], unavailable=check['explicit_unavailable'],
        rmsd_discrepancies=len(bad), short_mapping_discrepancies=len(bad_short),
        longer_mapping_discrepancies=len(bad-bad_short), short_alignments=short['alignments'],
        aligned_length_counts=lengths, maximum_rmsd_discrepancy=receipt['maximum_rmsd_rounding_error'],
        maximum_short_analytic_agreement_error=short['maximum_agreement_error'],
        scope='Full native-disposition/mapping/RMSD diagnostic and all short-case analytic checks completed. Original strict audit failure and discrepant native values remain quarantined. No complete geometry qualification, biological duplication effect or prediction accuracy is established.')
    with Path(plan['output']).open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
