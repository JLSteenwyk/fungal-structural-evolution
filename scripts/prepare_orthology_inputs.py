#!/usr/bin/env python3
"""Freeze a diverse 64-taxon reference core and all remaining taxa for assignment."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select(rows, qc):
    def quality(row):
        score = qc[row['taxon_id']]
        return (float(score['single_copy_percent']), float(score['complete_percent']),
                -float(score['duplicated_percent']), row['taxon_id'])
    fungi = [r for r in rows if r['study_role'] == 'ingroup']
    outgroups = [r for r in rows if r['study_role'] == 'outgroup']
    chosen, reasons = {}, {}
    # Broad taxonomic coverage is mandatory before quality ranks within a lineage.
    for role, candidates in [('ingroup', fungi), ('outgroup', outgroups)]:
        for group in sorted({r['lineage'].split(';')[0] for r in candidates}):
            row = max((r for r in candidates if r['lineage'].split(';')[0] == group), key=quality)
            chosen[row['taxon_id']] = row
            reasons[row['taxon_id']] = f'Core coverage of {role} group {group}'
    # Two extra close-holozoan representatives complement the eight outgroup categories.
    for group in ['Choanoflagellatea', 'Ichthyosporea']:
        row = max((r for r in outgroups if r['lineage'].split(';')[0] == group and r['taxon_id'] not in chosen), key=quality)
        chosen[row['taxon_id']] = row
        reasons[row['taxon_id']] = f'Additional core representation of {group}'
    if sum(r['study_role'] == 'outgroup' for r in chosen.values()) != 10:
        raise ValueError('Outgroup categories changed; review core design')

    def prefixes(row):
        parts = row['lineage'].split(';')
        return {tuple(parts[:i]) for i in range(2, len(parts) + 1)}

    while sum(r['study_role'] == 'ingroup' for r in chosen.values()) < 54:
        used = set().union(*(prefixes(r) for r in chosen.values() if r['study_role'] == 'ingroup'))
        counts = Counter(r['lineage'].split(';')[0] for r in chosen.values() if r['study_role'] == 'ingroup')
        def priority(row):
            novelty = sum(1 / len(prefix) for prefix in sorted(prefixes(row) - used))
            return (novelty, -counts[row['lineage'].split(';')[0]], *quality(row))
        row = max((r for r in fungi if r['taxon_id'] not in chosen), key=priority)
        chosen[row['taxon_id']] = row
        reasons[row['taxon_id']] = 'Additional core taxonomic-prefix diversity; QC tie-break'
    if len(chosen) != 64:
        raise ValueError('Unexpected core size')
    return chosen, reasons


def main():
    rows = read(ROOT / 'metadata/analysis_manifest.tsv')
    qc = {r['taxon_id']: r for r in read(ROOT / 'metadata/busco_eukaryota_qc.tsv')}
    source_receipt = ROOT / 'metadata/gene_representatives_receipt.json'
    receipt = json.loads(source_receipt.read_text())
    proteins = {r['taxon_id']: r for r in receipt['taxa']}
    if receipt['status'] != 'complete' or set(proteins) != {r['taxon_id'] for r in rows}:
        raise ValueError('Representative inputs do not cover the complete cohort')
    chosen, reasons = select(rows, qc)
    folder = ROOT / 'data/orthology_inputs/v1'
    if folder.exists():
        raise FileExistsError('Use an immutable new input version')
    for role in ['core', 'additional']:
        (folder / role).mkdir(parents=True)
    output = []
    for row in sorted(rows, key=lambda r: r['taxon_id']):
        name = row['taxon_id']
        source = ROOT / proteins[name]['path']
        if sha(source) != proteins[name]['sha256']:
            raise ValueError(f'Changed representative proteome: {name}')
        stage = 'core' if name in chosen else 'additional'
        target = folder / stage / f'{name}.faa'
        target.symlink_to(source.resolve())
        output.append({'taxon_id': name, 'species_name': row['species_name'], 'study_role': row['study_role'],
                       'lineage': row['lineage'], 'stage': stage, 'proteins': proteins[name]['selected_proteins'],
                       'source_sha256': proteins[name]['sha256'], 'input_path': str(target.relative_to(ROOT)),
                       'reason': reasons.get(name, 'Included in full-cohort species assignment and gene-tree analysis')})
    with (ROOT / 'metadata/orthology_input_manifest.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output)
    result = {'taxa': len(output), 'core_taxa': len(chosen), 'additional_taxa': len(output) - len(chosen),
              'core_proteins': sum(r['proteins'] for r in output if r['stage'] == 'core'),
              'additional_proteins': sum(r['proteins'] for r in output if r['stage'] == 'additional'),
              'core_fungal_groups': sorted({r['lineage'].split(';')[0] for r in chosen.values() if r['study_role'] == 'ingroup'}),
              'core_outgroup_groups': sorted({r['lineage'].split(';')[0] for r in chosen.values() if r['study_role'] == 'outgroup'}),
              'representative_receipt_sha256': sha(source_receipt),
              'input_manifest_sha256': sha(ROOT / 'metadata/orthology_input_manifest.tsv'),
              'selection_policy': 'Mandatory broad groups, then inverse-prefix-depth novelty; broad BUSCO single-copy recovery tie-break. Computational core, not a pilot.'}
    (ROOT / 'metadata/orthology_input_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
