#!/usr/bin/env python3
"""Partition recovered AlphaFold models into exact audited reuse and new ASA work."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from assess_pae_sensitivity import checked_receipt


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for name, digest in plan['pins'].items():
        assert sha(name) == digest, name
    full, previous, assessment, audit, output = [Path(plan[k]) for k in
        ['full', 'previous', 'assessment', 'audit', 'output']]
    assert not output.exists()
    for source in [full, previous, audit]:
        checked_receipt(source)
    ar = json.loads((audit / 'receipt.json').read_text())
    pr = json.loads((assessment / 'receipt.json').read_text())
    assert ar['status'] == 'passed_full_accessibility_snapshot'
    assert ar['snapshot_receipt_sha256'] == sha(previous / 'receipt.json')
    assert ar['assessment_receipt_sha256'] == sha(assessment / 'receipt.json')
    assert pr['status'] == 'complete_snapshot_predicted_accessibility'
    assert pr['config_sha256'] == sha(assessment / 'config.json')
    cfg = json.loads((assessment / 'config.json').read_text())
    assert cfg['snapshot_receipt_sha256'] == sha(previous / 'receipt.json')
    before = json.loads((previous / 'model_provenance.json').read_text())
    after = json.loads((full / 'model_provenance.json').read_text())
    old = {m['model_id']: m for m in before}
    new = {m['model_id']: m for m in after}
    assert len(old) == len(before) == ar['models_audited'] == pr['models']
    assert len(new) == len(after)
    assert sum(m['length'] for m in before) == ar['residues_audited'] == pr['residues']
    records = list(csv.DictReader((audit / 'audited_models.tsv').open(), delimiter='\t'))
    audited = {r['model_id']: r for r in records}
    assert len(audited) == len(records) and set(audited) == set(old) == set(pr['entry_receipts'])
    rows = []
    missing = []
    for sid in sorted(set(old) | set(new)):
        if sid not in new:
            rows.append(dict(model_id=sid, disposition='previous_only', length=old[sid]['length'],
                             metadata_differences='', entry_receipt_sha256='', table_sha256=''))
            continue
        if sid not in old:
            missing.append(new[sid])
            rows.append(dict(model_id=sid, disposition='new_accessibility_required', length=new[sid]['length'],
                             metadata_differences='', entry_receipt_sha256='', table_sha256=''))
            continue
        differences = {k for k in old[sid].keys() | new[sid].keys()
                       if old[sid].get(k) != new[sid].get(k)}
        # Only an explicitly documented extra provider identifier may differ.
        assert differences <= {'source_record_id'}, (sid, differences)
        if differences:
            assert 'source_record_id' not in old[sid]
            assert new[sid]['source_record_id'] == new[sid]['uniprot_accession']
        rp = assessment / (sid + '.receipt.json')
        rh = sha(rp)
        assert rh == pr['entry_receipts'][sid] == audited[sid]['entry_receipt_sha256']
        r = json.loads(rp.read_text())
        assert r['config_sha256'] == pr['config_sha256']
        assert r['model_sha256'] == new[sid]['sha256']
        assert r['sequence_sha256'] == new[sid]['sequence_sha256']
        assert r['residues'] == new[sid]['length']
        assert sha(assessment / (sid + '.residues.tsv.gz')) == r['table_sha256']
        rows.append(dict(model_id=sid, disposition='reuse_audited_accessibility', length=new[sid]['length'],
                         metadata_differences=','.join(sorted(differences)), entry_receipt_sha256=rh,
                         table_sha256=r['table_sha256']))
    output.mkdir(parents=True)
    (output / 'model_provenance.json').write_text(json.dumps(missing, indent=2) + '\n')
    with (output / 'model_dispositions.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    result = dict(status='complete_recovered_accessibility_gap_inventory', plan_sha256=sha(args.plan),
                  full_models=len(new), models=len(missing), residues=sum(m['length'] for m in missing),
                  reused_models=len(set(old) & set(new)), previous_only_models=len(set(old) - set(new)),
                  scope='Exact full-model partition; only additional source_record_id equal to UniProt accession allowed for reused models. Reused audited entry/table bytes checked. Coordinate bytes are bound by prior audit hashes, not reread here. New model_provenance is ASA input only, not a paired mapping.',
                  source_receipts={str(p / 'receipt.json'): sha(p / 'receipt.json') for p in [full, previous, assessment, audit]},
                  artifacts={p.name: sha(p) for p in output.iterdir()})
    (output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
