#!/usr/bin/env python3
"""Check every queue event, model, pair key and pair membership against frozen sources."""
import argparse
import csv
import json
import hashlib
from pathlib import Path
from estimate_expanded_duplication_pairs import sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ['queue', 'review', 'catalog', 'output']:
        ap.add_argument('--' + k, type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    rp = a.queue / 'receipt.json'
    r = json.loads(rp.read_text())
    assert r['status'] == 'complete_reviewed_duplication_model_pair_queue'
    sources = dict(r['source_hashes'])
    for root in [a.review, a.catalog]:
        assert str(root / 'receipt.json') in sources
    sources[str(rp)] = sha(rp)
    sources.update({str(a.queue / k): h for k, h in r['artifacts'].items()})
    for path, h in sources.items():
        assert sha(path) == h, path
    expected = {}
    for guide in ['profile', 'mafft']:
        path = a.review / (guide + '_candidate_tree_checks.tsv')
        assert str(path) in sources
        with path.open() as f:
            for row in csv.DictReader(f, delimiter='\t'):
                key = (guide, row['family'], row['gene_node'])
                assert key not in expected
                expected[key] = dict(row, guide=guide)
    models = {}
    with (a.queue / 'models.jsonl').open() as f:
        for line in f:
            row = json.loads(line)
            key = (row['model_id'], str(row['version']))
            assert key not in models
            models[key] = row
    needed_models = set()
    pair_keys = {}
    counts = dict(event_links=0, same_model_links=0, unresolved_tree_links=0, eligible_event_links=0)
    with (a.queue / 'event_model_pair_links.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            original = expected.pop((row['guide'], row['family'], row['gene_node']))
            assert all(row[k] == v for k, v in original.items())
            keys = sorted([(row['model_a'], row['version_a']), (row['model_b'], row['version_b'])])
            needed_models.update(keys)
            native = [(models[k]['model_id'], models[k]['version']) for k in keys]
            key = hashlib.sha256(json.dumps(native, separators=(',', ':')).encode()).hexdigest()
            assert row['pair_key'] == key
            counts['event_links'] += 1
            if row['tree_status'] != 'exact_reported_pair':
                status = 'unresolved_tree'
                counts['unresolved_tree_links'] += 1
            elif keys[0] == keys[1]:
                status = 'identical_model_no_alignment'
                counts['same_model_links'] += 1
            else:
                status = 'queued_distinct_models'
                counts['eligible_event_links'] += 1
                assert key not in pair_keys or pair_keys[key] == keys
                pair_keys[key] = keys
            assert row['comparison_status'] == status
    assert not expected and counts == r['counts'] and needed_models == set(models)
    seen = set()
    with (a.queue / 'model_pairs.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = row['pair_key']
            assert key not in seen
            seen.add(key)
            assert pair_keys[key] == [(row['model_a'], row['version_a']), (row['model_b'], row['version_b'])]
    assert seen == set(pair_keys)
    catalog_path = a.catalog / 'models.jsonl'
    assert str(catalog_path) in sources
    checked = set()
    with catalog_path.open() as f:
        for line in f:
            row = json.loads(line)
            key = (row['model_id'], str(row['version']))
            if key in models:
                assert key not in checked and models[key] == row
                checked.add(key)
    assert checked == set(models)
    active = {k for pair in pair_keys.values() for k in pair}
    assert len(models) == r['unique_models'] and len(active) == r['active_models']
    assert len(pair_keys) == r['unique_distinct_model_pairs'] and 2 * len(pair_keys) == r['directed_alignments_both_orders']
    assert sum(Path(models[k]['path']).stat().st_size for k in active) == r['active_coordinate_bytes']
    for path, h in sources.items():
        assert sha(path) == h, path
    result = dict(status='passed_full_duplication_model_pair_queue_export_readback',
                  source_sha256=sources, script_sha256=sha(__file__), counts=counts,
                  models_checked=len(models), distinct_pairs_checked=len(pair_keys),
                  scope='Every original reviewed event field, pair key/status/set, model membership and complete frozen catalog record checked. File sizes rechecked; raw coordinate bytes, protein-to-model source join, confidence qualification and alignment reuse eligibility are not independently reconstructed here.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_sha256', 'scope']}, indent=2))


if __name__ == '__main__':
    main()
