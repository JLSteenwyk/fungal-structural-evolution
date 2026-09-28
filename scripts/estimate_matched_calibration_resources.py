"""Inventory full-grid recorded fit costs; scenarios are not launch authorization."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from ancestral_chain_attempt import sha, write_json


def main():
    root = Path('results/structural_comparisons/full-matched-working-models-20260927-v1')
    out = Path('results/model_validation/matched-calibration-resources-20260928-v1')
    receipt = json.loads((root/'receipt.json').read_text())
    audit_path = Path('metadata/full_matched_working_model_output_audit_20260927.json')
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'passed_full_working_model_output_integrity_audit'
    assert audit['source_receipt_sha256'] == sha(root/'receipt.json')
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    rows = []; keys = set()
    for line in (root/'fit_manifest.jsonl').open():
        entry = json.loads(line); key = entry['fit_input_id'], entry['tree']
        assert key not in keys; keys.add(key)
        path = Path(entry['path']); assert sha(path) == entry['sha256']
        saved = json.loads(path.read_text()); p = saved['payload']
        assert (saved['fit_input_id'], saved['tree']) == key
        assert saved['plan_sha256'] == receipt['plan_sha256']
        seconds = p['wall_seconds']; assert math.isfinite(seconds) and seconds > 0
        rows.append(dict(fit_input_id=key[0], tree=key[1], records=p['records'],
                         recorded_fit_wall_seconds=seconds, serialized_bytes=path.stat().st_size,
                         status=p['status'], source_sha256=entry['sha256']))
    assert len(rows) == 144040
    frame = pd.DataFrame(rows)
    assert frame.fit_input_id.nunique() == 28808 and frame.tree.nunique() == 5
    total_seconds = float(frame.recorded_fit_wall_seconds.sum())
    total_bytes = int(frame.serialized_bytes.sum())
    scenarios = []
    for repeats in [199, 999, 1999]:
        for workers in [8, 16, 32]:
            scenarios.append(dict(replicates_per_unique_fit=repeats, workers=workers,
                total_refits=len(frame)*repeats,
                same_cost_serial_fit_hours=total_seconds*repeats/3600,
                ideal_parallel_days=total_seconds*repeats/3600/workers/24,
                same_format_fit_gib=total_bytes*repeats/2**30,
                coverage_mc_standard_error_at_095=math.sqrt(.95*.05/repeats),
                smallest_plus_one_tail_fraction=1/(repeats+1)))
    out.mkdir(parents=True, exist_ok=False)
    frame.to_parquet(out/'full_fit_costs.parquet', index=False)
    pd.testing.assert_frame_equal(pd.read_parquet(out/'full_fit_costs.parquet'), frame)
    pd.DataFrame(scenarios).to_csv(out/'resource_scenarios.tsv', sep='\t', index=False)
    result = dict(status='full_original_fit_cost_inventory_and_conditional_scenarios',
        unique_tree_fits=len(frame), unique_inputs=frame.fit_input_id.nunique(),
        recorded_fit_wall_seconds_sum=total_seconds, serialized_bytes_sum=total_bytes,
        fit_seconds_quantiles={str(q):float(np.quantile(frame.recorded_fit_wall_seconds,q)) for q in [0,.5,.9,.99,1]},
        source_receipt_sha256=sha(root/'receipt.json'), audit_sha256=sha(audit_path),
        script_sha256=sha(__file__), artifacts={p.name:sha(p) for p in out.iterdir()},
        scenarios=scenarios,
        scope='Full 144040 fit inventory. Recorded times are per-fit elapsed time, not CPU time. Linear scenarios assume identical refit cost and ideal parallelism; exclude simulation, input reconstruction, changed optimizer, retries, audits and contention. No production launch or calibrated inference.')
    write_json(out/'receipt.json', result)
    write_json(Path('metadata/matched_calibration_resources_20260928.json'), result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['artifacts','scenarios']}, indent=2))
    print(pd.DataFrame(scenarios).to_string(index=False))


if __name__ == '__main__':
    main()
