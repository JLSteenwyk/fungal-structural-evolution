#!/usr/bin/env python3
"""Integrate the complete ESMFold ASA cohorts with explicit relocation checks."""
import argparse
import json
import os
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def same_model(left,right):
    return left is not None and {k:v for k,v in left.items() if k!='path'}=={k:v for k,v in right.items() if k!='path'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', action='append', required=True, help='mapping=accessibility=audit')
    for name in ['snapshot', 'inventory', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--cohort-manifest',type=Path,required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use a new immutable union')
    mapping, inventory = checked_receipt(a.snapshot), checked_receipt(a.inventory)
    if inventory['status'] != 'complete_disjoint_prediction_inventory_union_with_batch_provenance' or mapping['source_inventory_sha256'] != sha(a.inventory / 'inventory.jsonl'):
        raise ValueError('Expected same-method combined model inventory')
    cohort_manifest=json.loads(a.cohort_manifest.read_text())
    if cohort_manifest['snapshot_receipt_sha256']!=sha(a.snapshot/'receipt.json'):raise ValueError('Wrong full snapshot')
    expected_sources={s['mapping']:s['mapping_receipt_sha256'] for s in cohort_manifest['cohorts']}
    if len(expected_sources)!=len(cohort_manifest['cohorts']):raise ValueError('Repeated manifest cohort')
    models = json.loads((a.snapshot / 'model_provenance.json').read_text())
    combined = {m['model_id']: m for m in models}
    entries = {}; summaries = []; sources = []; files = {}; seen = set(); settings = None
    for argument in a.cohort:
        mapping_name, asa_name, audit_name = argument.split('=')
        source, asa, audit = Path(mapping_name), Path(asa_name), Path(audit_name)
        mr, ar = checked_receipt(source), checked_receipt(audit)
        if mapping_name in seen or expected_sources.get(mapping_name) != sha(source / 'receipt.json'):
            raise ValueError('Unexpected or repeated source mapping')
        seen.add(mapping_name)
        receipt = json.loads((asa / 'receipt.json').read_text()); config = json.loads((asa / 'config.json').read_text())
        config_hash = sha(asa / 'config.json')
        if receipt['status'] != 'complete_snapshot_predicted_accessibility' or receipt['config_sha256'] != config_hash or config['snapshot_receipt_sha256'] != sha(source / 'receipt.json'):
            raise ValueError('Source ASA provenance differs')
        if ar['status'] != 'passed_full_accessibility_snapshot' or ar['assessment_receipt_sha256'] != sha(asa / 'receipt.json') or ar['snapshot_receipt_sha256'] != sha(source / 'receipt.json'):
            raise ValueError('Missing matching full raw source audit')
        method = {k: v for k, v in config.items() if k not in ['snapshot_receipt_sha256', 'models', 'residues', 'workers']}
        if settings is None:
            settings = method
        elif method != settings:
            raise ValueError('ASA method settings differ')
        provenance = json.loads((source / 'model_provenance.json').read_text())
        local = {m['model_id']: m for m in provenance}
        rows = read_table(audit / 'audited_models.tsv')
        if len(local) != len(provenance) or len(rows) != len(local) or {r['model_id'] for r in rows} != set(local) or set(receipt['entry_receipts']) != set(local):
            raise ValueError('Incomplete source model grid')
        if ar['models_audited'] != len(local) or ar['residues_audited'] != sum(m['length'] for m in local.values()):
            raise ValueError('Audit coverage differs')
        for row in rows:
            sid = row['model_id']; m = local[sid]
            if sid in entries or not same_model(combined.get(sid),m):
                raise ValueError('Overlapping or changed model provenance')
            rp = asa / (sid + '.receipt.json'); table = asa / (sid + '.residues.tsv.gz')
            digest = sha(rp); entry = json.loads(rp.read_text())
            if digest != row['entry_receipt_sha256'] or digest != receipt['entry_receipts'][sid] or entry['config_sha256'] != config_hash:
                raise ValueError('Changed entry or configuration')
            if entry['model_sha256'] != m['sha256'] or entry['sequence_sha256'] != m['sequence_sha256'] or int(row['residues']) != m['length']:
                raise ValueError('Model identity/length differs')
            if sha(table) != entry['table_sha256'] or sha(ROOT / m['path']) != m['sha256'] or sha(ROOT / combined[sid]['path']) != m['sha256']:
                raise ValueError('Changed audited residues/coordinates')
            entries[sid] = digest; summaries.append(row)
            files[rp.name] = rp; files[table.name] = table
        sources.append({'mapping': mapping_name, 'mapping_receipt_sha256': sha(source / 'receipt.json'),
                        'accessibility': asa_name, 'accessibility_receipt_sha256': sha(asa / 'receipt.json'),
                        'config_sha256': config_hash, 'workers': config['workers'], 'audit': audit_name, 'audit_receipt_sha256': sha(audit / 'receipt.json')})
    if seen != set(expected_sources) or set(entries) != set(combined) or len(entries) != mapping['distinct_models']:
        raise ValueError('Incomplete disjoint union')
    a.output.mkdir(parents=True)
    for name, source in files.items():
        (a.output / name).symlink_to(os.path.relpath(source.resolve(), a.output.resolve()))
    result = {'status': 'complete_disjoint_audited_accessibility_union', 'models': len(entries),
              'residues': sum(m['length'] for m in models), 'entry_receipts': entries,
              'mapping_receipt_sha256': sha(a.snapshot / 'receipt.json'), 'cohort_manifest_sha256':sha(a.cohort_manifest),'inventory_receipt_sha256':sha(a.inventory/'receipt.json'),'source_cohorts': sources,
              'shared_method_settings': settings, 'script_sha256': sha(Path(__file__)),
              'interpretation': 'Original entry receipts/residue tables are referenced through relative symlinks. Source raw audits and unchanged coordinate/table hashes establish union provenance; no new ASA calculation or independent raw-coordinate re-audit.'}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    audit = a.output / 'audit'; audit.mkdir()
    write_table(audit / 'audited_models.tsv', sorted(summaries, key=lambda r: r['model_id']))
    derived = {'status': 'passed_disjoint_accessibility_audit_union', 'models_audited': len(entries),
               'residues_audited': result['residues'], 'source_cohorts': sources,
               'assessment_receipt_sha256': sha(a.output / 'receipt.json'),
               'snapshot_receipt_sha256': sha(a.snapshot / 'receipt.json'),
               'script_sha256': sha(Path(__file__)), 'scope': result['interpretation'],
               'artifacts': {'audited_models.tsv': sha(audit / 'audited_models.tsv')}}
    (audit / 'receipt.json').write_text(json.dumps(derived, indent=2) + '\n')
    print(json.dumps(derived, indent=2), flush=True)


if __name__ == '__main__':
    main()
