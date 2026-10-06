#!/usr/bin/env python3
"""Independently read back a completed fresh MAFFT reconciliation recovery."""
import argparse
import csv
import json
from pathlib import Path

from audit_busco_gene_copies import ROOT, sha
from audit_hog_output_identities import audit as audit_hogs


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = read(args.plan)
    if sha(args.plan) != plan['plan_sha256']:
        raise ValueError('Readback plan changed')
    for name, digest in plan['pins'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Pinned source changed: ' + name)
    producer = ROOT / plan['producer_output']
    receipt = read(producer / 'receipt.json')
    if receipt['status'] != 'complete_native_mafft_reconciliation_pending_independent_output_readback':
        raise ValueError('Completed native producer receipt required')
    if receipt['guide'] != 'mafft' or receipt['statistics'] != {'species': 526, 'genes': 5815847}:
        raise ValueError('Unexpected producer scope')
    if receipt['input_manifest_sha256'] != plan['input_manifest_sha256']:
        raise ValueError('Producer input manifest differs from readback plan')
    result = ROOT / receipt['result']
    for name, digest in receipt['mandatory_artifacts'].items():
        if sha(result / name) != digest:
            raise ValueError('Mandatory native artifact changed: ' + name)
    with (result / 'Comparative_Genomics_Statistics/Statistics_Overall.tsv').open() as handle:
        stats = {row[0]: row[1] for row in csv.reader(handle, delimiter='\t') if len(row) >= 2}
    if int(stats['Number of species']) != 526 or int(stats['Number of genes']) != 5815847:
        raise ValueError('Independent statistics replay failed')

    output = ROOT / plan['output']
    if output.exists():
        raise FileExistsError('Readback output exists: ' + str(output))
    source = producer / 'Source/WorkingDirectory'
    audit_hogs(source, result, output / 'hog_identities')
    hog_receipt = read(output / 'hog_identities/receipt.json')
    if hog_receipt['status'] != 'passed_all_hog_identity_and_species_clade_checks':
        raise ValueError('HOG identity reader did not pass')
    for name, digest in receipt['mandatory_artifacts'].items():
        if sha(result / name) != digest:
            raise ValueError('Native artifact changed during independent readback: ' + name)
    result_receipt = {
        'status': 'passed_independent_mafft_reconciliation_recovery_output_readback',
        'producer_receipt_sha256': sha(producer / 'receipt.json'),
        'producer_output': plan['producer_output'],
        'guide': 'mafft',
        'statistics': {'species': int(stats['Number of species']), 'genes': int(stats['Number of genes'])},
        'hog_identity_receipt_sha256': sha(output / 'hog_identities/receipt.json'),
        'hog_node_tables': hog_receipt['node_tables'],
        'hog_rows': hog_receipt['hog_rows'],
        'gene_assignments_across_levels': hog_receipt['gene_assignments_across_levels'],
        'mandatory_artifacts': receipt['mandatory_artifacts'],
        'script_sha256': sha(__file__),
        'scope': 'Independent HOG identity, source-family, species-clade, table-universe, final-statistics and producer-artifact replay. It does not establish biological orthology, HOG ancestral completeness, duplication-event accuracy or guide robustness.',
    }
    (output / 'receipt.json').write_text(json.dumps(result_receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
