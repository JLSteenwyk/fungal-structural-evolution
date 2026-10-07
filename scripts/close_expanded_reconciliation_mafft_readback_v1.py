#!/usr/bin/env python3
"""Close a completed HOG audit after a final-receipt writer defect.

The original reader completed its all-node independent HOG audit and wrote the
HOG receipt, then failed only while hashing its own script path for the outer
receipt. This closure rechecks all protected sources, every HOG-table digest,
and the compact audit result before writing the missing outer receipt.
"""
import csv
import json
from pathlib import Path

from audit_busco_gene_copies import ROOT, sha

PLAN = ROOT / 'metadata/expanded_reconciliation_mafft_recovery_20261006_v2_readback_plan.json'
PRODUCER = ROOT / 'results/orthology/expanded-reconciliation-mafft-recovery-20261006-v2/receipt.json'
OUTPUT = ROOT / 'results/orthology/expanded-reconciliation-mafft-recovery-readback-20261006-v2'
HOG_OUTPUT = OUTPUT / 'hog_identities'
HOG_RECEIPT = HOG_OUTPUT / 'receipt.json'
SUMMARY = HOG_OUTPUT / 'node_summary.tsv'
CLOSURE = ROOT / 'metadata/expanded_reconciliation_mafft_recovery_20261006_v2_readback_closure.json'


def main():
    if (OUTPUT / 'receipt.json').exists() or CLOSURE.exists():
        raise FileExistsError('Immutable reconciliation readback closure already exists')
    plan = json.loads(PLAN.read_text())
    for name, digest in plan['pins'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Pinned source changed: ' + name)
    producer = json.loads(PRODUCER.read_text())
    if producer['status'] != 'complete_native_mafft_reconciliation_pending_independent_output_readback':
        raise ValueError('Completed native producer receipt required')
    result = ROOT / producer['result']
    for name, digest in producer['mandatory_artifacts'].items():
        if sha(result / name) != digest:
            raise ValueError('Native artifact changed: ' + name)
    with (result / 'Comparative_Genomics_Statistics/Statistics_Overall.tsv').open() as handle:
        stats = {row[0]: row[1] for row in csv.reader(handle, delimiter='\t') if len(row) >= 2}
    if (int(stats['Number of species']), int(stats['Number of genes'])) != (526, 5815847):
        raise ValueError('Native statistics replay failed')
    hog = json.loads(HOG_RECEIPT.read_text())
    if hog['status'] != 'passed_all_hog_identity_and_species_clade_checks':
        raise ValueError('All-node independent HOG audit did not pass')
    if (hog['source_proteins'], hog['source_families'], hog['node_tables']) != (5815847, 658522, 525):
        raise ValueError('Unexpected independent HOG audit scope')
    if sha(SUMMARY) != hog['artifacts']['node_summary.tsv']:
        raise ValueError('HOG summary changed')
    with SUMMARY.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 525 or {row['node'] for row in rows} != {'N' + str(i) for i in range(525)}:
        raise ValueError('HOG summary node universe differs')
    folder = result / 'Phylogenetic_Hierarchical_Orthogroups'
    for row in rows:
        path = folder / (row['node'] + '.tsv')
        if not path.is_file() or sha(path) != row['sha256']:
            raise ValueError('HOG table changed after independent audit: ' + row['node'])
    for path, digest in hog['input_hashes'].items():
        if sha(Path(path)) != digest:
            raise ValueError('HOG audit source changed: ' + path)
    receipt = {
        'status': 'passed_independent_mafft_reconciliation_recovery_output_readback',
        'producer_receipt_sha256': sha(PRODUCER),
        'producer_output': str(Path(producer['result'])),
        'guide': 'mafft',
        'statistics': {'species': int(stats['Number of species']), 'genes': int(stats['Number of genes'])},
        'hog_identity_receipt_sha256': sha(HOG_RECEIPT),
        'hog_node_tables': hog['node_tables'],
        'hog_rows': hog['hog_rows'],
        'gene_assignments_across_levels': hog['gene_assignments_across_levels'],
        'mandatory_artifacts': producer['mandatory_artifacts'],
        'closure_script_sha256': sha(Path(__file__)),
        'writer_defect': "The original reader's all-node audit completed, then its final receipt writer passed a string rather than Path to sha(). This closure independently rehashes every HOG table and source/summary artifact before adding the missing outer receipt.",
        'scope': 'Independent HOG identity, source-family, species-clade, table-universe, final-statistics and producer-artifact replay. It does not establish biological orthology, HOG ancestral completeness, duplication-event accuracy or guide robustness.'
    }
    closure = {
        'status': 'passed_post_audit_readback_receipt_closure_after_writer_defect',
        'readback_receipt_sha256': None,
        'hog_receipt_sha256': receipt['hog_identity_receipt_sha256'],
        'producer_receipt_sha256': receipt['producer_receipt_sha256'],
        'all_hog_tables_rehashed': len(rows),
        'source_hashes': {str(p.relative_to(ROOT)): sha(p) for p in (PLAN, PRODUCER, HOG_RECEIPT, SUMMARY, Path(__file__))},
        'scope': 'Closure of a post-audit writer defect; records the independent all-node audit and rehash closure. No producer or HOG-audit artifact is altered.'
    }
    (OUTPUT / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    closure['readback_receipt_sha256'] = sha(OUTPUT / 'receipt.json')
    CLOSURE.write_text(json.dumps(closure, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'hog_node_tables': receipt['hog_node_tables'], 'hog_rows': receipt['hog_rows']}, indent=2))


if __name__ == '__main__':
    main()
