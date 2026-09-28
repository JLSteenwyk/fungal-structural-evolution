#!/usr/bin/env python3
"""Account for all 324 diagnostic configurations, preserving failed attempts."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path('results/ancestral')
PASS = 'independent_tip_tree_and_candidate_readback_passed'


def main():
    pins = {}

    def sha(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def read(path):
        path = Path(path)
        pins[str(path)] = sha(path)
        return json.loads(path.read_text())

    def audit(folder):
        directory = ROOT/folder
        receipt = read(directory/'receipt.json')
        for name, digest in receipt['artifacts'].items():
            assert sha(directory/name) == digest
        for name, digest in receipt['pins'].items():
            assert sha(name) == digest
        return read(directory/'dispositions.json')

    base = audit('historian-capacity-audit-20260927-final-v1')
    assert len(base)==324 and len({r['job_id'] for r in base})==324
    retries = {}
    for method, folder in [('mafft','historian-memory-retry-final-readback-20260927-v1'),
                           ('famsa','historian-famsa-memory-retry-final-readback-20260927-v1')]:
        producer = ROOT/('historian-memory-retries-20260927-v1' if method=='mafft' else 'historian-famsa-memory-retries-20260927-v1')
        for row in audit(folder):
            assert row['status']==PASS and row['job_id'] not in retries
            retries[row['job_id']] = (row,producer)
    assert len(retries)==4
    outrows = []
    for original in base:
        jid = original['job_id']
        original_path = ROOT/'historian-capacity-20260927-v1'/jid/'receipt.json'
        attempt = read(original_path)
        selected = original
        selected_path = original_path
        if original['status'] != PASS:
            assert original['status']=='recorded_unsuccessful_execution'
            selected,producer = retries[jid]
            selected_path = producer/jid/'receipt.json'
            retry = read(selected_path)
            assert retry['job'] == attempt['job']
            assert retry['exit_code']==0 and attempt['exit_code']!=0
            for name,digest in retry['artifacts'].items():
                assert sha(selected_path.parent/name)==digest
        else:
            assert jid not in retries
        outrows.append(dict(job_id=jid,input_id=original['input_id'],original_status=original['status'],
            selected_status=selected['status'],used_memory_retry=jid in retries,
            original_receipt=str(original_path),selected_receipt=str(selected_path),
            selected_elapsed_seconds=selected['elapsed_seconds'],
            selected_peak_sampled_rss_bytes=selected['peak_sampled_rss_bytes']))
    assert sum(r['used_memory_retry'] for r in outrows)==4
    comparison = read('metadata/historian_largest_family_comparison_completed_20260927.json')
    assert sha(comparison['completed_receipt_path'])==comparison['completed_receipt_sha256']
    summary = read(ROOT/'historian-capacity-sensitivity-summary-20260927-final-v1/receipt.json')
    assert summary['dispositions']=={PASS:320,'recorded_unsuccessful_execution':4}
    output=ROOT/'historian-integrated-capacity-20260927-v1'
    output.mkdir(exist_ok=False)
    with (output/'configuration_dispositions.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,list(outrows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(outrows)
    receipt=dict(status='all_324_configurations_have_independently_checked_diagnostic_outputs',
        configurations=324,original_passes=320,original_failed_attempts_preserved=4,checked_memory_retries=4,
        original_grid_dependent_comparisons=summary['dependent_node_comparisons'],
        original_grid_changed_comparisons=summary['changed_node_comparisons'],
        largest_family_nonroot_sensitivity=comparison['nonroot_summary'],
        pins=pins,script_sha256=sha(__file__),
        artifacts={'configuration_dispositions.tsv':sha(output/'configuration_dispositions.tsv')},
        scope='Execution and output integrity complete for this fixed-parameter diagnostic grid. Original failures retained alongside successful exact-input retries. Original-grid zero-edit comparisons exclude failed attempts; retry sensitivity must be reported separately. No posterior convergence, model adequacy, historical root identification or qualified ancestral ensemble is established.')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    completion=dict(receipt,completed_receipt_path=str(output/'receipt.json'),completed_receipt_sha256=sha(output/'receipt.json'))
    Path('metadata/historian_integrated_capacity_completed_20260927.json').write_text(json.dumps(completion,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':
    main()
