"""Frozen full triad/input sources for sequence-only correspondence controls."""
import json
from pathlib import Path

from full_triad_fit_sources import load_sources as load_mapping
from reference_measurement_union_sources import bind, verify


def load_original(plan, plan_path):
    original_path = Path(plan['geometry_source_plan'])
    original = json.loads(original_path.read_text())
    ready, inputs, mapping, root, bindings = load_mapping(original, original_path)
    for path, digest in plan['pins'].items():
        bind(bindings, path, digest)
    bind(bindings, plan_path)
    mp = json.loads(Path(original['mapping_plan']).read_text())
    tp = json.loads(Path(mp['triad_plan']).read_text())
    path = Path(tp['output']) / 'ordered_model_triads.jsonl'
    catalog = [json.loads(line) for line in path.open()]
    assert len(catalog) == 31235 and len(ready) == 27056
    assert [t for t in catalog if t['source_design_ready_links']] == ready
    assert mapping['target_contexts'] == 283409
    assert mapping['reference_tie_records'] == 214461
    assert mapping['duplicate_reference_links'] == 428922
    # Full sequences determine correspondence before confidence filtering.
    for t in ready:
        for model in t['models']:
            r = inputs[(*model, 'full')]
            assert r['status'] == 'ready'
            assert len(r['sequence']) == r['original_length']
            assert r['original_positions'] == list(range(1, r['original_length'] + 1))
    verify(bindings)
    return catalog, ready, inputs, mapping, bindings


def load_catalog(plan, plan_path):
    bindings = dict(plan['pins'])
    bind(bindings, plan_path)
    closed = json.loads(Path(plan['catalog_completion']).read_text())
    assert closed['status'] == 'complete_verified_full_triad_sequence_catalog'
    assert closed['ordered_model_triads'] == 31235
    assert closed['source_ready_triads'] == 27056
    assert closed['exact_process_journals_checked'] == 2
    bind(bindings, plan['catalog_completion'])
    bind(bindings, closed['full_hash_archive'], closed['full_hash_archive_sha256'])
    archive = json.loads(Path(closed['full_hash_archive']).read_text())
    assert len(archive['services']) == 2
    for path, digest in archive['source_hashes'].items():
        bind(bindings, path, digest)
    root = Path(closed['producer_receipt']).parent
    sets = [json.loads(line) for line in (root / 'sequence_sets.jsonl').open()]
    assert len(sets) == closed['unique_sequence_model_sets']
    assert 12 * len(sets) == closed['native_alignment_states']
    verify(bindings)
    return sets, closed, bindings
