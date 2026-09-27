"""Bind full discrepancy accounting and analytic short-mapping checks."""
import json
from pathlib import Path
import subprocess
from screen_duplication_domain_alignment_coverage import sha


def main():
    units = ['fungal-primary-alignment-rmsd-diagnostic-20260927.service',
             'fungal-primary-rmsd-diagnostic-completion-20260927.service',
             'fungal-primary-short-geometry-readback-20260927.service']
    states = {}
    for unit in units:
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), (unit, state)
        states[unit] = state
    root = Path('results/structural_comparisons/primary-alignment-rmsd-diagnostic-20260927-v1')
    rp = root/'receipt.json'
    cp = Path('metadata/primary_alignment_rmsd_diagnostic_completed_20260927.json')
    sp = Path('metadata/primary_alignment_short_geometry_readback_20260927.json')
    receipt, check, short = [json.loads(p.read_text()) for p in [rp, cp, sp]]
    assert receipt['status'] == 'complete_primary_alignment_rmsd_diagnostic_not_scientific_acceptance'
    assert receipt['scientific_eligibility'] is False
    assert check['status'] == 'verified_completed_primary_rmsd_diagnostic_not_scientific_acceptance'
    assert short['status'] == 'passed_all_short_primary_alignment_analytic_rmsd_checks'
    assert check['receipt_sha256'] == short['source_receipt_sha256'] == sha(rp)
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    assert check['numeric_rows'] == receipt['numerically_checked_alignments']
    lengths = check['aligned_length_counts']
    assert short['alignments'] == sum(lengths.get(str(n), 0) for n in [1, 2])
    assert short['one_residue'] == lengths.get('1', 0)
    assert short['two_residue'] == lengths.get('2', 0)
    bad = {(r['pair_key'], r['mask'], r['order']) for r in receipt['rmsd_discrepancies']}
    bad_short = {(r['pair_key'], r['mask'], r['order']) for r in short['rows'] if r['rmsd_status'] == 'outside_printed_rounding'}
    assert bad_short <= bad and len(bad_short) == short['outside_printed_rounding']
    result = dict(status='complete_verified_primary_rmsd_diagnostic_and_short_geometry',
        scientific_eligibility=False, terminal_states=states, source_hashes={str(p):sha(p) for p in [rp, cp, sp]},
        script_sha256=sha(__file__), dispositions=receipt['directed_dispositions'],
        numeric_alignments=check['numeric_rows'], unavailable=check['explicit_unavailable'],
        rmsd_discrepancies=len(bad), short_mapping_discrepancies=len(bad_short),
        longer_mapping_discrepancies=len(bad-bad_short), short_alignments=short['alignments'],
        aligned_length_counts=lengths, maximum_rmsd_discrepancy=receipt['maximum_rmsd_rounding_error'],
        maximum_short_analytic_agreement_error=short['maximum_agreement_error'],
        scope='Full native-disposition/mapping/RMSD diagnostic and all short-case analytic checks completed. Original strict audit failure and discrepant native values remain quarantined. No complete geometry qualification, biological duplication effect or prediction accuracy is established.')
    Path('metadata/primary_rmsd_diagnostic_and_short_completed_20260927.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
