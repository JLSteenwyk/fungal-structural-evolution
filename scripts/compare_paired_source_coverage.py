#!/usr/bin/env python3
"""Compare coverage contributions without pooling predictor-specific characters."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha
from prepare_paired_phylogenetic_inputs import write_table


def read_source(path):
    receipt = checked_receipt(path)
    if receipt['status'] != 'complete_paired_phylogenetic_input_preparation':
        raise ValueError('Completed paired inputs required')
    with (path/'marker_summary.tsv').open() as f:
        markers = {r['marker']: r for r in csv.DictReader(f, delimiter='\t')}
    with (path/'taxon_coverage.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    grid = {}
    for row in rows:
        key = row['taxon_id'], row['marker']
        if key in grid or key[1] not in markers:
            raise ValueError('Repeated or unknown coverage key')
        grid[key] = row['taxon_eligible'] == 'True' and markers[key[1]]['status'] == 'ready_for_inference'
    if len(grid) != receipt['taxon_marker_rows']:
        raise ValueError('Coverage grid count differs')
    return receipt, grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['reference', 'local', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new immutable comparison output')
    rr, reference = read_source(args.reference)
    lr, local = read_source(args.local)
    if rr['source_receipts']['matrix']['sha256'] != lr['source_receipts']['matrix']['sha256'] or set(reference) != set(local):
        raise ValueError('Sources must share the exact alignment and coverage grid')
    if rr['mask'] != lr['mask'] or rr['eligibility'] != lr['eligibility']:
        raise ValueError('Sources use different eligibility definitions')
    manifest = ROOT/'metadata/analysis_manifest.tsv'
    with manifest.open() as f:
        taxa = {r['taxon_id']: r for r in csv.DictReader(f, delimiter='\t')}
    matrix = Path(rr['source_receipts']['matrix']['path'])
    if sha(matrix/'receipt.json') != rr['source_receipts']['matrix']['sha256']:
        raise ValueError('Changed matrix receipt')
    if json.loads((matrix/'receipt.json').read_text())['manifest_sha256'] != sha(manifest):
        raise ValueError('Changed manifest')
    marker_ids = {m for t,m in reference}
    if set(reference) != {(t,m) for t in taxa for m in marker_ids}:
        raise ValueError('Incomplete manifest coverage grid')
    rows = []
    for t, taxon in taxa.items():
        counts = Counter()
        for m in marker_ids:
            a,b = reference[t,m],local[t,m]
            counts['both' if a and b else 'reference_only' if a else 'local_only' if b else 'neither'] += 1
        rows.append({'taxon_id':t, 'species_name':taxon['species_name'], 'study_role':taxon['study_role'],
                     'lineage_group':taxon['lineage'].split(';')[0],
                     **{k:counts[k] for k in ['reference_only','local_only','both','neither']},
                     'reference_usable':counts['reference_only']+counts['both'],
                     'local_usable':counts['local_only']+counts['both'],
                     'union_usable':len(marker_ids)-counts['neither']})
    args.output.mkdir(parents=True)
    write_table(args.output/'taxon_source_contributions.tsv', rows)
    receipt = {'status':'complete_predictor_source_coverage_comparison',
        'reference_receipt_sha256':sha(args.reference/'receipt.json'),
        'local_receipt_sha256':sha(args.local/'receipt.json'),
        'reference_path':str(args.reference), 'local_path':str(args.local),
        'manifest_sha256':sha(manifest), 'script_sha256':sha(Path(__file__)),
        'taxa':len(taxa), 'markers':len(marker_ids),
        'taxa_with_reference':sum(r['reference_usable']>0 for r in rows),
        'taxa_with_local':sum(r['local_usable']>0 for r in rows),
        'taxa_with_either':sum(r['union_usable']>0 for r in rows),
        'newly_represented_taxa':sum(r['reference_usable']==0 and r['local_usable']>0 for r in rows),
        'taxon_marker_counts':{k:sum(r[k] for r in rows) for k in ['reference_only','local_only','both','neither']},
        'interpretation':'Coverage union only; no pooled alignments or branch estimates. Each usable combination independently meets its source-specific four-taxon family eligibility. This can undercount combinations enabled by a future pooled analysis. Unequal predictor calibration, orthology and tree support remain unresolved.',
        'artifacts':{'taxon_source_contributions.tsv':sha(args.output/'taxon_source_contributions.tsv')}}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
