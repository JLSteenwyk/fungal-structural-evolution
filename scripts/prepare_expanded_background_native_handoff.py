#!/usr/bin/env python3
"""Normalize the complete current input union without selecting favorable results."""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path
from expanded_background_input_sources import load_sources, MASKS, PREFERENCE
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def digest(value): return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run(plan_path):
    plan_path = Path(plan_path); plan = json.loads(plan_path.read_text())
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib'] * 2 ** 30
    active, pairs, candidates, specifications, bindings = load_sources(plan, plan_path)
    owners = {}; model_counts = {}; models = {}
    for label, source in specifications.items():
        seen = set()
        with source['model_path'].open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version']; assert key not in seen; seen.add(key)
                fingerprint = [row[k] for k in ['sha256', 'sequence_sha256', 'length']]
                if key in owners: assert owners[key]['fingerprint'] == fingerprint
                state = owners.setdefault(key, dict(fingerprint=fingerprint, labels=[])); state['labels'].append(label)
                if key in active:
                    assert fingerprint == [active[key][k] for k in ['sha256', 'sequence_sha256', 'length']]
                    models[key] = active[key]
        assert len(seen) == source['receipt']['models']; model_counts[label] = len(seen)
    assert set(models) == set(active)
    overlapping = {key for key, value in owners.items() if len(value['labels']) > 1}
    semantics = {}; selected = {}; counts = Counter(); collection_states = {}; actual_pdb = set(); written_bytes = 0
    out = Path(plan['output']); out.mkdir(exist_ok=False, parents=True)
    with gzip.open(out / 'overlapping_input_states.jsonl.gz', 'xt', compresslevel=1) as overlap_output:
        for label, source in specifications.items():
            seen = set(); source_counts = Counter()
            with source['manifest'].open() as handle:
                for ordinal, line in enumerate(handle, 1):
                    row = json.loads(line); key = row['model_id'], row['version'], row['mask']
                    assert key not in seen and key[2] in MASKS and label in owners[key[:2]]['labels']; seen.add(key)
                    assert row['source_sha256'] == owners[key[:2]]['fingerprint'][0]
                    assert row['status'] in ['ready', 'too_few_retained_residues', 'source_rejected']; source_counts[key[2] + ':' + row['status']] += 1
                    normalized = {k: v for k, v in row.items() if k not in ['path', 'coordinate_shard']}; signature = digest(normalized)
                    if key[:2] in overlapping:
                        assert key not in semantics or semantics[key] == signature; semantics[key] = signature
                        overlap_output.write(json.dumps(dict(collection=label, source_manifest=str(source['manifest']), source_row_number=ordinal,
                                                             source_row_sha256=digest(row), semantic_sha256=signature, input=row), separators=(',', ':')) + '\n')
                    if key[:2] not in active and key[:2] not in overlapping: continue
                    if row['status'] == 'ready':
                        assert len(row['sequence']) == row['retained_residues'] == len(row['original_positions'])
                        assert len(set(row['original_positions'])) == len(row['original_positions']) and all(1 <= p <= row['original_length'] for p in row['original_positions'])
                        if row['path'] not in actual_pdb:
                            bind(bindings, row['path'], row['sha256']); assert sha(row['path']) == row['sha256']; actual_pdb.add(row['path']); written_bytes += Path(row['path']).stat().st_size
                    if key[:2] in active:
                        origin = dict(collection=label, source_manifest=str(source['manifest']), source_row_number=ordinal, source_row_sha256=digest(row),
                                      path=row.get('path'), sha256=row.get('sha256'), coordinate_shard=row.get('coordinate_shard'))
                        if key in selected:
                            assert selected[key]['semantic_sha256'] == signature; selected[key]['collection_sources'].append(origin)
                        else:
                            keep = {k: row[k] for k in ['model_id', 'version', 'mask', 'status', 'source_sha256', 'sequence', 'original_length', 'retained_residues', 'reason', 'path', 'sha256'] if k in row}
                            selected[key] = dict(**keep, semantic_sha256=signature, collection_sources=[origin])
            assert len(seen) == source['receipt']['input_dispositions'] == 2 * model_counts[label] and dict(source_counts) == source['receipt']['counts']
            assert {key[:2] for key in seen} == {key for key in owners if label in owners[key]['labels']}
            collection_states[label] = len(seen)
            print('Full input collection verified', label, len(seen), flush=True)
    assert set(selected) == {(*key, mask) for key in active for mask in MASKS}
    assert set(semantics) == {(*key, mask) for key in overlapping for mask in MASKS}
    raw_bytes = 0
    with gzip.open(out / 'active_models.jsonl.gz', 'xt', compresslevel=1) as handle:
        for key, row in sorted(models.items()):
            bind(bindings, row['path'], row['sha256']); assert sha(row['path']) == row['sha256']; raw_bytes += Path(row['path']).stat().st_size
            handle.write(json.dumps(row, separators=(',', ':')) + '\n')
    with gzip.open(out / 'active_inputs.jsonl.gz', 'xt', compresslevel=1) as handle:
        for key, row in sorted(selected.items()):
            assert row['collection_sources'][0]['collection'] == min(owners[key[:2]]['labels'], key=PREFERENCE.index)
            counts[key[2] + ':' + row['status']] += 1; handle.write(json.dumps(row, separators=(',', ':')) + '\n')
    partition = []; seen_pairs = set(); native_counts = Counter(); new = pending = 0
    for pair in sorted(pairs, key=lambda row: row['pair_key']):
        ends = [(pair['model_a'], int(pair['version_a'])), (pair['model_b'], int(pair['version_b']))]; key = pair['pair_key']; prior = candidates[key]
        assert ends[0] != ends[1] and ends == sorted(ends) and key not in seen_pairs; seen_pairs.add(key)
        assert hashlib.sha256(json.dumps(ends, separators=(',', ':')).encode()).hexdigest() == key
        assert ends == [(prior['model_a'], int(prior['version_a'])), (prior['model_b'], int(prior['version_b']))]
        assert all((*end, mask) in selected for end in ends for mask in MASKS)
        matching = json.loads(prior['matching_old_sources'])
        if matching:
            assert prior['disposition'] == 'matching_catalog_sources_pending_input_and_result_checks'; disposition = 'pending_actual_input_result_and_numeric_reuse_checks'; pending += 1
        else:
            assert prior['disposition'] in ['new_pair_requires_alignment', 'changed_catalog_sources_require_alignment']; disposition = 'native_measurement_required'; new += 1
            for mask in MASKS:
                ready = all(selected[(*end, mask)]['status'] == 'ready' for end in ends)
                native_counts[mask + ':' + ('ready' if ready else 'input_unavailable')] += 2
        partition.append(dict(**pair, **{k: prior[k] for k in ['new_sources', 'matching_old_sources', 'changed_old_sources']}, catalog_disposition=prior['disposition'], measurement_disposition=disposition))
    assert seen_pairs == set(candidates) and new == plan['expected']['new_pairs'] and pending == plan['expected']['pending_catalog_reuse_pairs']
    with (out / 'full_background_work_partition.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(partition[0]), delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(partition)
    summary = dict(active_models=len(active), active_model_mask_states=len(selected), full_pairs=len(pairs), new_pairs=new, pending_catalog_reuse_pairs=pending,
                   collection_model_counts=model_counts, collection_input_states=collection_states, overlapping_models=len(overlapping), overlapping_model_mask_states=len(semantics),
                   active_input_counts=dict(counts), directed_native_input_counts=dict(native_counts), actual_ready_PDB_files=len(actual_pdb), actual_raw_coordinate_models=len(models), raw_coordinate_bytes=raw_bytes, ready_PDB_bytes=written_bytes)
    verify(bindings)
    names = ['active_models.jsonl.gz', 'active_inputs.jsonl.gz', 'overlapping_input_states.jsonl.gz', 'full_background_work_partition.tsv']
    result = dict(status='complete_full_expanded_background_native_input_handoff_pending_independent_readback', plan_sha256=sha(plan_path), **summary,
                  input_collection_preference=PREFERENCE, source_hashes=bindings, artifacts={name: sha(out / name) for name in names}, scientific_eligibility=False, scope=plan['scope'])
    with (out / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', required=True, type=Path); run(parser.parse_args().plan)
