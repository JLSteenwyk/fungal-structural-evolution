#!/usr/bin/env python3
"""Verify every full-reference original length and both complete frozen input manifests."""
import argparse
import json
from pathlib import Path
from reference_order_coverage_sources import load_inventory_inputs
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    args = p.parse_args(); plan = json.loads(args.plan.read_text()); out = Path(plan['output']); assert not out.exists()
    native = json.loads(Path(plan['native_plan']).read_text()); bindings = {str(args.plan): sha(args.plan), **plan['pins'], **native['pins']}
    bind(bindings, plan['native_plan']); verify(bindings)
    pairs, models, inputs = load_inventory_inputs(native, plan, bindings)
    verify(bindings)
    result = dict(status='passed_full_reference_original_length_and_written_input_preflight',
                  plan_sha256=sha(args.plan), native_plan_sha256=sha(plan['native_plan']), full_pairs=len(pairs), models=len(models),
                  reference_input_dispositions=len(inputs), all_input_dispositions=plan['expected']['all_input_dispositions'],
                  source_hashes=bindings, scientific_eligibility=False,
                  scope='Full inventory/readback and both entire frozen written-input manifests verified. Every reference model/version/raw source hash/original full length, full sequence hash, retained count/position/sequence length and both masks checked; duplicate/missing states rejected. Native closed union remains required before order/coverage production. No alignment or geometry rerun and no scientific eligibility acceptance.')
    with out.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
