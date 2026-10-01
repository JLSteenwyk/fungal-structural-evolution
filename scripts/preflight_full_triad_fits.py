#!/usr/bin/env python3
"""Measure full mapping work and storage before same-residue geometric fits."""
import argparse
import json
from collections import Counter
from pathlib import Path
from full_triad_fit_sources import load_sources, iterate_maps, DEFINITIONS, triple_sha
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text())
    triads, inputs, r, root, bindings = load_sources(plan, args.plan)
    counts = Counter(); needed = set(); states = occurrences = unique = unique_residues = 0
    previous = None; seen = set()
    for item in iterate_maps(root, triads, inputs):
        if item['triad_id'] != previous:
            previous = item['triad_id']; seen = set()
        for definition, field in DEFINITIONS:
            status = item['common_core_fit_input_status'][definition]; counts[item['mask'] + ':' + definition + ':' + status] += 1
            triples = item[field]; occurrences += len(triples)
            if status != 'pending_common_coordinate_geometry': continue
            key = item['mask'], triple_sha(triples)
            if key not in seen:
                seen.add(key); unique += 1; unique_residues += len(triples)
            needed.update((*model, item['mask']) for model in item['models'])
        states += 1
        if states % 50000 == 0: print('Full common-fit work preflight', states, '/', r['mask_order_states'], flush=True)
    assert states == r['mask_order_states'] and occurrences == r['common_reference_residue_occurrences'] + r['cycle_consistent_residue_occurrences']
    pdb_bytes = sum(Path(inputs[key]['path']).stat().st_size for key in needed)
    assert all(inputs[key]['status'] == 'ready' for key in needed)
    verify(bindings)
    result = dict(status='passed_full_triad_geometric_fit_work_preflight', plan_sha256=sha(args.plan), mapping_states=states,
                  fit_rows=states * 2, counts=dict(counts), triple_occurrences=occurrences,
                  unique_eligible_triple_sets_within_physical_triad_and_mask=unique,
                  unique_eligible_residue_occurrences=unique_residues,
                  proper_pair_fits_per_implementation=3 * unique, needed_pdb_inputs=len(needed), needed_pdb_bytes=pdb_bytes,
                  source_hashes=bindings, scientific_eligibility=False,
                  scope='Exhaustive full closed mapping/input grid, both definitions/masks/all8orders. Counts actual eligible coordinate work and within-triad/mask identical-triple reuse; no sampled geometry, preselection or pilot. RAM/disk/runtime estimates in immutable resource plan remain estimates, not observed runtime or scientific qualification.')
    with args.output.open('x') as f: f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
