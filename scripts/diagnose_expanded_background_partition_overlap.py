#!/usr/bin/env python3
"""Exhaust model-partition overlaps and compare all overlapping written inputs."""
import json
from itertools import combinations
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    root = Path('results/structural_comparisons'); bindings = {__file__: sha(__file__)}
    specifications = dict(primary=('duplication-model-pair-queue-20260928-v1', 'models.jsonl'),
                          old_reference_additional=('duplication-reference-comparison-inventory-20260926-v1', 'additional_models.jsonl'),
                          expanded_background_additional=('background-measurement-inventory-20260930-v1', 'additional_models.jsonl'))
    catalogs = {}
    for label, (folder, filename) in specifications.items():
        rp = root / folder / 'receipt.json'; receipt = json.loads(rp.read_text()); bind(bindings, rp)
        path = rp.parent / filename; bind(bindings, path, receipt['artifacts'][filename]); catalog = {}
        with path.open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version']; assert key not in catalog; catalog[key] = row['sha256']
        catalogs[label] = catalog
    intersections = []; overlaps = set()
    for left, right in combinations(catalogs, 2):
        keys = catalogs[left].keys() & catalogs[right].keys()
        assert all(catalogs[left][key] == catalogs[right][key] for key in keys)
        intersections.append(dict(partitions=[left, right], overlapping_models=len(keys), different_source_hashes=0)); overlaps.update(keys)
    assert len(overlaps) == 6
    inputs = {}; input_folders = dict(primary='duplication-alignment-inputs-20260928-v1', old_reference_additional='duplication-reference-alignment-inputs-20260926-v1')
    for label, folder in input_folders.items():
        rp = root / folder / 'receipt.json'; receipt = json.loads(rp.read_text()); bind(bindings, rp)
        manifest = rp.parent / 'inputs.jsonl'; bind(bindings, manifest, receipt['artifacts']['inputs.jsonl']); selected = {}
        with manifest.open() as handle:
            for line in handle:
                row = json.loads(line); key = row['model_id'], row['version'], row['mask']
                if key[:2] in overlaps: assert key not in selected; selected[key] = row
        assert set(selected) == {(*key, mask) for key in overlaps for mask in ['full', 'plddt70']}; inputs[label] = selected
    comparisons = []
    for key in sorted(inputs['primary']):
        left, right = [inputs[label][key] for label in input_folders]
        # Physical fields must be identical; preserve collection-specific paths
        # and coordinate-shard provenance separately rather than treating them as
        # physical differences or silently dropping them from a future union.
        ignored = {'path', 'coordinate_shard'}
        assert {k: v for k, v in left.items() if k not in ignored} == {k: v for k, v in right.items() if k not in ignored}
        assert left['source_sha256'] == catalogs['primary'][key[:2]] == catalogs['old_reference_additional'][key[:2]]
        assert left['status'] == right['status'] == 'ready' and left['sha256'] == right['sha256']
        for row in [left, right]: bind(bindings, row['path'], row['sha256'])
        comparisons.append(dict(model_mask=list(key), status='ready', source_sha256=left['source_sha256'], written_pdb_sha256=left['sha256'],
                                input_paths=[left['path'], right['path']], coordinate_shards=[left['coordinate_shard'], right['coordinate_shard']],
                                original_length=left['original_length'], retained_residues=left['retained_residues'], semantic_fields_equal=True))
    verify(bindings)
    result = dict(status='passed_full_partition_overlap_and_all_overlapping_written_input_diagnostic', model_partition_sizes={k: len(v) for k, v in catalogs.items()},
                  intersections=intersections, overlapping_model_mask_states=len(comparisons), written_pdb_files_checked=2 * len(comparisons), comparisons=comparisons,
                  source_hashes=bindings, alignment_reuse_qualified=False, complete_input_union_qualified=False, scientific_eligibility=False,
                  scope='Exhaustive intersection of all three model partitions, then every full/p70 state for all six overlapping models. Exact semantic input fields, source model hashes and all24actual written PDB file hashes agree; distinct collection paths/shards retained. This diagnoses why the old strict nonoverlap handoff cannot be used unchanged for expanded backgrounds. Does not qualify the full active-model handoff, any old alignment reuse/native geometry, matched effect or scientific inference.')
    path = Path('metadata/expanded_background_partition_overlap_diagnostic_20261001.json')
    with path.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['source_hashes', 'comparisons']}, indent=2))


if __name__ == '__main__': main()
