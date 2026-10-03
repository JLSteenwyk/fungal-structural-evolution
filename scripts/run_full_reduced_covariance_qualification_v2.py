#!/usr/bin/env python3
"""Full retained-basis qualification after an explicit parallel source closure."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile

from ancestral_chain_attempt import sha
from full_covariance_qualification_sources import jsonl, LINK_EXTRA
from full_exact_covariance_sources import runtime_caps
from full_reduced_covariance_sources_v2 import load, NEW_LINK_EXTRA, SUMMARY_FIELDS
from reduced_covariance_basis import SCHEMA, reduce_record, readback_record, identity
from reference_measurement_union_sources import verify, bind

PRODUCER_STATUS = 'complete_full_exact_retained_uniform_covariance_qualification_pending_readback_v2'
READER_STATUS = 'passed_full_exact_retained_uniform_covariance_qualification_independent_readback_v2'
ARTIFACTS = ['stage_plan.json', 'design_covariance_audits.jsonl.gz', 'setting_audit_links.tsv.gz']


def database(path):
    db = sqlite3.connect(path)
    db.execute('CREATE TABLE audits(old_id TEXT PRIMARY KEY,aid TEXT UNIQUE,cohort TEXT,did TEXT,mode TEXT,tree TEXT,status TEXT,old_status TEXT,UNIQUE(did,mode,tree))')
    return db


def add(db, old, row):
    db.execute('INSERT INTO audits VALUES(?,?,?,?,?,?,?,?)',
        (old['audit_id'], row['audit_id'], row['cohort_id'], row['design_id'],
        row['loading_mode'], row['tree'], row['disposition'], old['disposition']))


def finish(source, db, counts, linkcounts, basis, audits, links, settings):
    old = source['old']; ncohorts = db.execute('SELECT COUNT(DISTINCT cohort) FROM audits').fetchone()[0]
    ndesigns = db.execute('SELECT COUNT(DISTINCT did) FROM audits').fetchone()[0]
    assert db.execute('SELECT COUNT(*) FROM audits').fetchone()[0] == audits == old['audit_rows']
    assert [r[0] for r in db.execute('SELECT DISTINCT cohort FROM audits ORDER BY cohort')] == source['cohort_ids']
    assert ncohorts == old['unique_cohorts'] and ndesigns == old['unique_designs']
    assert links == old['setting_audit_links'] and settings == old['model_setting_rows']
    return dict(logical_cases=old['logical_cases'], model_setting_rows=settings,
        unique_cohorts=ncohorts, unique_designs=ndesigns, audit_rows=audits, setting_audit_links=links,
        audit_status_counts=dict(counts), link_status_counts=dict(linkcounts), trees=old['trees'],
        loading_modes=old['loading_modes'], retained_basis_audit_counts=dict(basis))


def run(path, reader=False):
    path = Path(path); plan = json.loads(path.read_text())
    runtime_caps(); source, bindings = load(plan, path); original_bindings = dict(bindings)
    root = Path(plan['output']); receipt_path = root / 'receipt.json'
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    marker = dict(plan_sha256=sha(path), schema=SCHEMA, source_contract=source['contract'])
    if reader:
        receipt = json.loads(receipt_path.read_text()); assert receipt['status'] == PRODUCER_STATUS
        assert receipt['plan_sha256'] == sha(path) and receipt['source_contract'] == source['contract']
        assert receipt['source_hashes'] == original_bindings and receipt['scientific_eligibility'] is False
        assert set(receipt['artifacts']) == set(ARTIFACTS)
        bind(bindings, receipt_path)
        for name, expected in receipt['artifacts'].items():
            assert sha(root / name) == expected; bind(bindings, root / name, expected)
        assert json.loads((root / 'stage_plan.json').read_text()) == marker
        assert not (root / 'readback.json').exists()
    else:
        root.mkdir(exist_ok=False)
        with (root / 'stage_plan.json').open('x') as f: f.write(json.dumps(marker, indent=2) + '\n')
    counts = Counter(); linkcounts = Counter(); basis = Counter(); audits = links = settings = 0
    old_records = jsonl(source['old_root'] / 'design_covariance_audits.jsonl.gz')
    with tempfile.TemporaryDirectory(prefix='retained-audit-index-', dir=root) as directory:
        db = database(str(Path(directory) / 'audits.sqlite'))
        if reader:
            exported = jsonl(root / 'design_covariance_audits.jsonl.gz')
            for old in old_records:
                row = next(exported)
                cert = source['certificates'][old['cohort_id'], old['loading_mode']]
                readback_record(row, old, cert, source['contract'])
                assert row['tree'] in source['old']['trees']
                add(db, old, row); counts[row['disposition']] += 1
                basis[str(len(row['retained_kernel_names']))] += 1; audits += 1
                if audits % 30000 == 0:
                    db.commit(); print('independent_retained_covariance_audits', audits, '/', source['old']['audit_rows'], flush=True)
            assert next(exported, None) is None
        else:
            with gzip.open(root / 'design_covariance_audits.jsonl.gz', 'wt') as f:
                for old in old_records:
                    cert = source['certificates'][old['cohort_id'], old['loading_mode']]
                    row = reduce_record(old, cert, source['contract'])
                    assert row['tree'] in source['old']['trees']
                    add(db, old, row); counts[row['disposition']] += 1
                    basis[str(len(row['retained_kernel_names']))] += 1; audits += 1
                    f.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
                    if audits % 30000 == 0:
                        db.commit(); print('retained_covariance_audits', audits, '/', source['old']['audit_rows'], flush=True)
        db.commit()
        with gzip.open(source['old_root'] / 'setting_audit_links.tsv.gz', 'rt') as f:
            old_links = csv.DictReader(f, delimiter='\t'); fields = source['settings_fields']
            assert old_links.fieldnames == fields + LINK_EXTRA
            if reader:
                target = gzip.open(root / 'setting_audit_links.tsv.gz', 'rt')
                exported_links = csv.DictReader(target, delimiter='\t')
                assert exported_links.fieldnames == fields + NEW_LINK_EXTRA
            else:
                target = gzip.open(root / 'setting_audit_links.tsv.gz', 'wt')
                writer = csv.DictWriter(target, fields + NEW_LINK_EXTRA, delimiter='\t', lineterminator='\n'); writer.writeheader()
            with target:
                for old in old_links:
                    hit = db.execute('SELECT aid,did,mode,tree,status,old_status,cohort FROM audits WHERE old_id=?', (old['audit_id'],)).fetchone()
                    assert hit is not None
                    aid, did, mode, tree, status, old_status, cohort = hit
                    assert (did, mode, tree, old_status, cohort) == (old['design_id'], old['loading_mode'], old['tree'], old['covariance_disposition'], old['cohort_id'])
                    assert aid == identity(source['contract'], old['audit_id'])
                    combined = old['disposition'] if old['disposition'] != 'ready_for_working_covariance_fit' else status
                    if reader:
                        row = next(exported_links)
                        # Independent field checks retain every original setting,
                        # outcome, fit-input identity and non-ready disposition.
                        assert {k: row[k] for k in fields} == {k: old[k] for k in fields}
                        assert [row[k] for k in NEW_LINK_EXTRA] == [mode, tree, aid, status, combined, old['audit_id'], old_status]
                    else:
                        writer.writerow({**{k: old[k] for k in fields}, 'loading_mode': mode, 'tree': tree,
                            'audit_id': aid, 'covariance_disposition': status, 'combined_disposition': combined,
                            'source_audit_id': old['audit_id'], 'source_covariance_disposition': old_status})
                    links += 1; linkcounts[combined] += 1
                    if links % (len(source['old']['trees']) * len(source['old']['loading_modes'])) == 0: settings += 1
                if reader: assert next(exported_links, None) is None
        summary = finish(source, db, counts, linkcounts, basis, audits, links, settings); db.close()
    verify(bindings)
    result = dict(status=READER_STATUS if reader else PRODUCER_STATUS, plan_sha256=sha(path),
        source_contract=source['contract'], **summary, source_hashes=bindings,
        scientific_eligibility=False, production_fitting_launched=False, scope=plan['scope'])
    if reader:
        assert all(receipt[k] == summary[k] for k in SUMMARY_FIELDS)
        result['producer_receipt_sha256'] = sha(receipt_path)
        output = root / 'readback.json'
    else:
        result['artifacts'] = {name: sha(root / name) for name in ARTIFACTS}; output = receipt_path
    with output.open('x') as f: f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--reader', action='store_true'); a = p.parse_args(); run(a.plan, a.reader)
