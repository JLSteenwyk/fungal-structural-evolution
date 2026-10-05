#!/usr/bin/env python3
"""Observe exact original full domain-reader and atlas-union controllers."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from record_process_covariance_pipeline_checkpoint import observe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); assert not args.output.exists()
    handles = []; states = {}; pins = {str(Path(__file__)): sha(__file__)}
    for prefix, memory in [('completed_afdb_domain_registry_readback', 16), ('full_prediction_atlas_union', 32)]:
        launch = Path(f'metadata/{prefix}_launch_20261005_v1.json')
        record = json.loads(launch.read_text()); plan = Path(record['plan'])
        expected = {'cpu.max':'200000 100000','memory.max':str(memory*2**30),'memory.swap.max':'0'}
        handles.append(observe(str(launch), expected))
        assert sha(plan) == record['plan_sha256']
        pins[str(launch)] = sha(launch); pins[str(plan)] = record['plan_sha256']
        root = Path(json.loads(plan.read_text())['output'])
        state = root/'state.json'
        states[prefix] = dict(output_root_present=root.exists(),
            own_final_receipt_present=Path(f'metadata/{prefix}_20261005_v1.json').exists())
        if state.exists(): states[prefix]['checkpoint'] = json.loads(state.read_text())
    result = dict(status='verified_original_full_atlas_update_runtime',
        checked_utc=datetime.now(timezone.utc).isoformat(),original_handles=handles,states=states,
        source_hashes=pins,full_taxa=526,full_representative_proteins=5815847,
        scientific_eligibility=False,gpu=False,new_predictions=0,all_eight_aims_incomplete=True,
        scope='Exact original PID/create/command/cgroup or invocation-linked terminal journal '
              'and immutable launch/plan hashes. Live output counts are checkpoints only; this '
              'does not qualify full domain arithmetic, atlas union, confidence or biological aims.')
    with args.output.open('x') as handle: json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','original_handles']},indent=2))


if __name__ == '__main__': main()
