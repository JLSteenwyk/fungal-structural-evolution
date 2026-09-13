#!/usr/bin/env python3
"""Prepare expanded paired evolutionary inputs after verified confidence integration."""
import argparse
import fcntl
import json
import subprocess
import sys
import time
from pathlib import Path
from advance_expanded_structural_confidence import process_identity
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--confidence-pid', type=int, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    control = ROOT / 'results/structural_alphabet/expanded-confidence-control-v1'
    qualified = ROOT / 'results/structural_alphabet/audited-gdm-expanded-v1'
    mapping = ROOT / 'results/structural_markers/gdm-expanded-v1'
    matrix = ROOT / 'results/phylogeny/profile-matrix-50-v1'
    baseline_metadata = ROOT / 'metadata/paired_phylogenetic_input_receipt.json'
    baseline_receipt = json.loads(baseline_metadata.read_text())
    # Locate the actual baseline output by its versioned receipt, never by newest mtime.
    matches = [x.parent for x in (ROOT/'results').glob('**/receipt.json')
               if x.parent.name.startswith('paired') and sha(x) == sha(baseline_metadata)]
    if len(matches) != 1:
        raise ValueError('Expected one exact baseline paired input snapshot')
    baseline = matches[0]
    checked_receipt(baseline)
    inputs = ROOT / 'results/phylogeny/paired-inputs-gdm-expanded-v1'
    a.output.mkdir(parents=True, exist_ok=True)
    lock = (a.output/'.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    current = process_identity(a.confidence_pid)
    cp = a.output/'config.json'
    if current and 'scripts/advance_expanded_structural_confidence.py' not in current['command']:
        raise ValueError('Unexpected dependency process')
    dependency = json.loads(cp.read_text())['dependency_identity'] if cp.exists() else current
    scripts = ['run_expanded_paired_inputs.py', 'prepare_paired_phylogenetic_inputs.py',
               'advance_expanded_structural_confidence.py', 'assess_pae_sensitivity.py',
               'compare_marker_structures.py', 'audit_busco_gene_copies.py']
    config = {'dependency_pid': a.confidence_pid, 'dependency_identity': dependency,
              'confidence_config_sha256': sha(control/'config.json'),
              'mapping_receipt_sha256': sha(mapping/'receipt.json'),
              'matrix_receipt_sha256': sha(matrix/'receipt.json'),
              'baseline_path': str(baseline.relative_to(ROOT)),
              'baseline_receipt_sha256': sha(baseline_metadata),
              'scripts': {n: sha(ROOT/'scripts'/n) for n in scripts},
              'inputs_output': str(inputs.relative_to(ROOT))}
    if cp.exists() and json.loads(cp.read_text()) != config:
        raise ValueError('Changed expanded input configuration')
    cp.write_text(json.dumps(config, indent=2)+'\n')
    last = 0
    while current:
        if current != dependency:
            raise ValueError('Dependency PID reused')
        if time.monotonic()-last >= 60:
            print('Waiting for confidence integration PID', a.confidence_pid, flush=True)
            last = time.monotonic()
        time.sleep(10)
        current = process_identity(a.confidence_pid)
    cr = json.loads((control/'receipt.json').read_text())
    if cr['status'] != 'complete_expanded_mapping_bound_structural_confidence' or cr['config_sha256'] != config['confidence_config_sha256'] or cr['qualified_receipt_sha256'] != sha(qualified/'receipt.json'):
        raise ValueError('Complete verified confidence integration required')
    for name, digest in config['scripts'].items():
        if sha(ROOT/'scripts'/name) != digest:
            raise ValueError('Changed producer '+name)
    for directory, digest in [(mapping, config['mapping_receipt_sha256']), (matrix, config['matrix_receipt_sha256'])]:
        if sha(directory/'receipt.json') != digest:
            raise ValueError('Changed paired source')
    if not (inputs/'receipt.json').exists():
        command = [sys.executable, str(ROOT/'scripts/prepare_paired_phylogenetic_inputs.py'),
                   '--encodings', str(qualified), '--snapshot', str(mapping),
                   '--matrix', str(matrix), '--output', str(inputs)]
        print('Starting expanded paired alignments', flush=True)
        with (a.output/'preparation.log').open('a') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True, cwd=ROOT)
    r = checked_receipt(inputs)
    if r['status'] != 'complete_paired_phylogenetic_input_preparation' or r['source_receipts']['encodings']['sha256'] != sha(qualified/'receipt.json'):
        raise ValueError('Completed expanded inputs differ')
    if r['mask'] != baseline_receipt['mask'] or r['eligibility'] != baseline_receipt['eligibility'] or r['source_receipts']['matrix']['sha256'] != baseline_receipt['source_receipts']['matrix']['sha256']:
        raise ValueError('Coverage snapshots have incomparable alignment/selection rules')
    old = {x['marker']: x for x in read_table(baseline/'marker_summary.tsv')}
    new = {x['marker']: x for x in read_table(inputs/'marker_summary.tsv')}
    if set(old) != set(new):
        raise ValueError('Marker universe differs')
    rows = [{'marker': k, 'previous_eligible_taxa': old[k]['eligible_taxa'],
             'expanded_eligible_taxa': new[k]['eligible_taxa'],
             'eligible_taxon_count_change': int(new[k]['eligible_taxa'])-int(old[k]['eligible_taxa']),
             'previous_retained_columns': old[k]['retained_columns'],
             'expanded_retained_columns': new[k]['retained_columns'],
             'previous_status': old[k]['status'], 'expanded_status': new[k]['status']}
            for k in sorted(new)]
    write_table(a.output/'marker_coverage_change.tsv', rows)
    result = {'status': 'complete_expanded_paired_inputs_and_coverage_comparison',
              'config_sha256': sha(cp), 'confidence_receipt_sha256': sha(control/'receipt.json'),
              'input_receipt_sha256': sha(inputs/'receipt.json'),
              'taxa_audited': r['taxa_audited'], 'markers_audited': r['markers_audited'],
              'previous_ready_markers': baseline_receipt['ready_markers'], 'expanded_ready_markers': r['ready_markers'],
              'artifacts': {'marker_coverage_change.tsv': sha(a.output/'marker_coverage_change.tsv')},
              'interpretation': 'Identical matrix and confidence/eligibility rules across acquisition snapshots; source availability and chosen models can change. Counts are not evidence of evolutionary acceleration. Gene-copy, alignment and prediction-source sensitivity remain necessary.'}
    (a.output/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
