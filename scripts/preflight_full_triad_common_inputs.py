#!/usr/bin/env python3
"""Validate full original/masked input grid and bound potential common-map output."""
import argparse
import json
from pathlib import Path
from full_triad_common_sources import load_design_inputs
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    triads, inputs, design, counts, bindings = load_design_inputs(plan); bind(bindings, args.plan)
    upper_triples = 0; max_length = 0; lengths = []; input_counts = {}
    for triad in triads:
        for mask in ['full', 'plddt70']:
            rows = [inputs[(*model, mask)] for model in triad['models']]
            upper_triples += 16 * min(r['retained_residues'] for r in rows)  # 8 orders, both definitions
            max_length = max(max_length, *(r['original_length'] for r in rows))
    for row in inputs.values():
        key = row['mask'] + ':' + row['status']; input_counts[key] = input_counts.get(key, 0) + 1
        if row['mask'] == 'full': lengths.append(row['original_length'])
    lengths.sort(); verify(bindings)
    result = dict(status='passed_full_triad_original_and_masked_input_preflight', plan_sha256=sha(args.plan),
                  target_contexts=design['target_contexts'], reference_tie_records=design['reference_tie_records'], duplicate_reference_links=design['duplicate_reference_links'],
                  unique_ordered_model_triads=design['unique_ordered_model_triads'], correspondence_work_triads=len(triads), potential_correspondence_states=16 * len(triads),
                  needed_models=len(inputs) // 2, all_input_dispositions=plan['expected']['all_input_dispositions'], input_source_counts=counts, needed_input_counts=input_counts,
                  maximum_original_length=max_length, minimum_original_length=lengths[0], median_original_length=lengths[len(lengths) // 2],
                  upper_bound_triple_occurrences_both_definitions=upper_triples,
                  upper_bound_uncompressed_triple_bytes=upper_triples * (3 * len(str(max_length)) + 5),
                  source_hashes=bindings, scientific_eligibility=False,
                  scope='Actual complete closed triad design and both full written-input manifests checked. Required full sequences/source hashes/original lengths/nonconsecutive positions/masked subsequences/statuses and full input universes verified. Potential triple/JSON-size upper bounds use actual retained lengths, both masks/all8orders and both common/cycle definitions; they exclude metadata/headers/compression overhead. Not native mapping, geometric fitting, coverage-cohort or biological qualification.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
