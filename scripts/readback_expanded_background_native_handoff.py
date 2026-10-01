#!/usr/bin/env python3
"""Independently rebuild the full current background input union with SQL."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path
from expanded_background_input_sources import load_sources, MASKS, PREFERENCE, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path, output = Path(plan_path), Path(output); assert not output.exists(); plan = json.loads(plan_path.read_text())
    active, pairs, candidates, specifications, expected_bindings = load_sources(plan, plan_path)
    root = Path(plan['output']); rp = root / 'receipt.json'; r = json.loads(rp.read_text())
    assert r['status'] == 'complete_full_expanded_background_native_input_handoff_pending_independent_readback' and r['plan_sha256'] == sha(plan_path)
    assert r['input_collection_preference'] == PREFERENCE and r['scientific_eligibility'] is False and r['scope'] == plan['scope']
    assert set(r['artifacts']) == {'active_models.jsonl.gz', 'active_inputs.jsonl.gz', 'overlapping_input_states.jsonl.gz', 'full_background_work_partition.tsv'}
    bindings = dict(r['source_hashes']); bind(bindings, rp)
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE models(m TEXT,v INT,c TEXT,raw TEXT,seq TEXT,n INT,PRIMARY KEY(m,v,c))')
    db.execute('CREATE TABLE inputs(m TEXT,v INT,mask TEXT,rank INT,signature TEXT,payload TEXT,origin TEXT,status TEXT,path TEXT,sha TEXT,PRIMARY KEY(m,v,mask,rank))')
    db.execute('CREATE TABLE chosen(m TEXT,v INT,mask TEXT,status TEXT,PRIMARY KEY(m,v,mask))')
    model_counts = {}; collection_states = {}
    for label, source in specifications.items():
        count = 0; batch = []
        with source['model_path'].open() as handle:
            for line in handle:
                row = json.loads(line); batch.append((row['model_id'], row['version'], label, row['sha256'], row['sequence_sha256'], row['length'])); count += 1
                if len(batch) == 10000: db.executemany('INSERT INTO models VALUES(?,?,?,?,?,?)', batch); batch.clear()
        if batch: db.executemany('INSERT INTO models VALUES(?,?,?,?,?,?)', batch)
        assert count == source['receipt']['models']; model_counts[label] = count
    assert db.execute('SELECT COUNT(*) FROM (SELECT m,v FROM models GROUP BY m,v HAVING COUNT(DISTINCT raw)!=1 OR COUNT(DISTINCT seq)!=1 OR COUNT(DISTINCT n)!=1)').fetchone()[0] == 0
    overlapping = {(m, v) for m, v in db.execute('SELECT m,v FROM models GROUP BY m,v HAVING COUNT(*)>1')}
    overlap_groups = set(); actual_pdb = set(); written_bytes = 0
    with gzip.open(root / 'overlapping_input_states.jsonl.gz', 'rt') as exported_overlaps:
        for label, source in specifications.items():
            count = 0; seen = set(); statuses = Counter(); batch = []
            with source['manifest'].open() as handle:
                for ordinal, line in enumerate(handle, 1):
                    row = json.loads(line); key = row['model_id'], row['version'], row['mask']; assert key not in seen and key[2] in MASKS; seen.add(key); count += 1
                    catalog = db.execute('SELECT raw FROM models WHERE m=? AND v=? AND c=?', (*key[:2], label)).fetchone(); assert catalog and catalog[0] == row['source_sha256']
                    assert row['status'] in ['ready', 'too_few_retained_residues', 'source_rejected']; statuses[key[2] + ':' + row['status']] += 1
                    projection = {k: row[k] for k in row if k not in {'path', 'coordinate_shard'}}
                    signature = hashlib.sha256(json.dumps(projection, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                    row_hash = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                    if key[:2] in overlapping:
                        line = next(exported_overlaps, None); assert line is not None
                        actual = json.loads(line); expected = dict(collection=label, source_manifest=str(source['manifest']), source_row_number=ordinal, source_row_sha256=row_hash, semantic_sha256=signature, input=row)
                        assert actual == expected; overlap_groups.add(key)
                    if key[:2] not in active and key[:2] not in overlapping: continue
                    if row['status'] == 'ready':
                        assert row['retained_residues'] == len(row['sequence']) == len(row['original_positions'])
                        assert len(set(row['original_positions'])) == row['retained_residues'] and all(1 <= p <= row['original_length'] for p in row['original_positions'])
                        if row['path'] not in actual_pdb:
                            bind(expected_bindings, row['path'], row['sha256']); assert bindings[row['path']] == row['sha256']; actual_pdb.add(row['path']); written_bytes += Path(row['path']).stat().st_size
                    compact = {k: row[k] for k in ['model_id', 'version', 'mask', 'status', 'source_sha256', 'sequence', 'original_length', 'retained_residues', 'reason', 'path', 'sha256'] if k in row}
                    origin = dict(collection=label, source_manifest=str(source['manifest']), source_row_number=ordinal, source_row_sha256=row_hash,
                                  path=row.get('path'), sha256=row.get('sha256'), coordinate_shard=row.get('coordinate_shard'))
                    batch.append((*key, PREFERENCE.index(label), signature, json.dumps(compact), json.dumps(origin), row['status'], row.get('path'), row.get('sha256')))
                    if len(batch) == 10000: db.executemany('INSERT INTO inputs VALUES(?,?,?,?,?,?,?,?,?,?)', batch); batch.clear()
            if batch: db.executemany('INSERT INTO inputs VALUES(?,?,?,?,?,?,?,?,?,?)', batch)
            assert count == source['receipt']['input_dispositions'] == 2 * model_counts[label] and dict(statuses) == source['receipt']['counts']
            assert {key[:2] for key in seen} == {(m, v) for m, v in db.execute('SELECT m,v FROM models WHERE c=?', (label,))}; collection_states[label] = count
            print('Independent full input collection', label, count, flush=True)
        assert next(exported_overlaps, None) is None
    assert db.execute('SELECT COUNT(*) FROM (SELECT m,v,mask FROM inputs GROUP BY m,v,mask HAVING COUNT(DISTINCT signature)!=1)').fetchone()[0] == 0
    assert overlap_groups == {(*key, mask) for key in overlapping for mask in MASKS}
    counts = Counter(); states = 0; raw_bytes = 0
    with gzip.open(root / 'active_models.jsonl.gz', 'rt') as models, gzip.open(root / 'active_inputs.jsonl.gz', 'rt') as inputs:
        for key, original in sorted(active.items()):
            actual = json.loads(next(models)); assert actual == original
            catalog = db.execute('SELECT raw,seq,n FROM models WHERE m=? AND v=?', key).fetchall(); assert catalog
            assert all((original['sha256'], original['sequence_sha256'], original['length']) == row for row in catalog)
            bind(expected_bindings, original['path'], original['sha256']); assert bindings[original['path']] == original['sha256']; raw_bytes += Path(original['path']).stat().st_size
            for mask in MASKS:
                rows = list(db.execute('SELECT signature,payload,origin FROM inputs WHERE m=? AND v=? AND mask=? ORDER BY rank', (*key, mask))); assert rows
                expected = json.loads(rows[0][1]); expected.update(semantic_sha256=rows[0][0], collection_sources=[json.loads(row[2]) for row in rows])
                actual = json.loads(next(inputs)); assert actual == expected and all(row[0] == rows[0][0] for row in rows)
                db.execute('INSERT INTO chosen VALUES(?,?,?,?)', (*key, mask, actual['status'])); counts[mask + ':' + actual['status']] += 1; states += 1
        assert next(models, None) is next(inputs, None) is None
    assert states == 2 * len(active) and r['source_hashes'] == expected_bindings
    pending = new = 0; native_counts = Counter(); seen_pairs = set()
    with (root / 'full_background_work_partition.tsv').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for source in sorted(pairs, key=lambda row: row['pair_key']):
            actual = next(reader, None); assert actual is not None; key = source['pair_key']; prior = candidates[key]
            ends = [(source['model_a'], int(source['version_a'])), (source['model_b'], int(source['version_b']))]
            assert ends == sorted(ends) and ends[0] != ends[1] and key not in seen_pairs; seen_pairs.add(key)
            assert hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest() == key
            assert ends == [(prior['model_a'], int(prior['version_a'])), (prior['model_b'], int(prior['version_b']))]
            matching = json.loads(prior['matching_old_sources'])
            if matching:
                assert prior['disposition'] == 'matching_catalog_sources_pending_input_and_result_checks'; disposition = 'pending_actual_input_result_and_numeric_reuse_checks'; pending += 1
            else:
                assert prior['disposition'] in ['new_pair_requires_alignment', 'changed_catalog_sources_require_alignment']; disposition = 'native_measurement_required'; new += 1
            expected = {**source, **{k: prior[k] for k in ['new_sources', 'matching_old_sources', 'changed_old_sources']}, 'catalog_disposition': prior['disposition'], 'measurement_disposition': disposition}
            assert actual == expected
            for mask in MASKS:
                statuses = [db.execute('SELECT status FROM chosen WHERE m=? AND v=? AND mask=?', (*end, mask)).fetchone() for end in ends]; assert all(statuses)
                if not matching: native_counts[mask + ':' + ('ready' if all(row[0] == 'ready' for row in statuses) else 'input_unavailable')] += 2
        assert next(reader, None) is None
    assert seen_pairs == set(candidates) and new == plan['expected']['new_pairs'] and pending == plan['expected']['pending_catalog_reuse_pairs']
    summary = dict(active_models=len(active), active_model_mask_states=states, full_pairs=len(pairs), new_pairs=new, pending_catalog_reuse_pairs=pending,
                   collection_model_counts=model_counts, collection_input_states=collection_states, overlapping_models=len(overlapping), overlapping_model_mask_states=len(overlap_groups),
                   active_input_counts=dict(counts), directed_native_input_counts=dict(native_counts), actual_ready_PDB_files=len(actual_pdb), actual_raw_coordinate_models=len(active),
                   raw_coordinate_bytes=raw_bytes, ready_PDB_bytes=written_bytes)
    assert set(summary) == set(SUMMARY_FIELDS) and all(r[k] == value for k, value in summary.items()); verify(bindings)
    result = dict(status='passed_full_expanded_background_native_input_handoff_sql_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=sha(rp), checker_sha256=sha(__file__),
                  **summary, source_hashes=bindings, scientific_eligibility=False,
                  scope='Entire four current/legacy catalog/input universes independently indexed with SQL; every global overlap requires identical full physical semantic hashes. Every normalized active model/mask/origin/source preference/nullable state, raw/PDB binding, complete background work partition and new-native input disposition independently reconstructed. Source/proof I/O shared, no producer normalization or partition code imported. Does not qualify old-result reuse, native alignment geometry or biological matched effects.')
    with output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.plan, args.output)
