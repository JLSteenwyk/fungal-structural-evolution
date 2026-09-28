"""Read back every serialized interval task against the selected export."""
import json
from pathlib import Path
import pandas as pd
from ancestral_chain_attempt import sha, write_json

root = Path('results/model_validation/selected-matched-interval-grid-inputs-20260928-v1')
r = json.loads((root/'receipt.json').read_text())
for name, digest in r['artifacts'].items():
    assert sha(root/name) == digest
for name, digest in r['source_bindings'].items():
    assert sha(name) == digest
selected = pd.read_parquet('results/structural_comparisons/refined-working-model-grid-export-20260928-v2/unique_fits.parquet')
expected = {(row.fit_input_id, row.tree): row.selected_source_sha256 for row in selected.itertuples()}
seen = set()
inputs = set()
for line in (root/'tasks.jsonl').open():
    task = json.loads(line)
    identifier = task['fit_input_id']
    assert identifier not in inputs
    inputs.add(identifier)
    assert task['cache_entry']['fit_input_id'] == identifier
    assert set(task['selected_sources']) == set(r['factors'])
    for tree, source in task['selected_sources'].items():
        key = (identifier, tree)
        assert key not in seen and expected[key] == source['sha256']
        seen.add(key)
assert seen == set(expected) and len(inputs) == 28808 and len(seen) == 144040
result = dict(status='passed_all_serialized_selected_interval_task_bindings',
    inputs=len(inputs), fits=len(seen), preparation_receipt_sha256=sha(root/'receipt.json'),
    script_sha256=sha(__file__), scope='Full task/source mapping readback; interval execution and cache replay still separate prerequisites.')
write_json(Path('metadata/selected_matched_interval_grid_readback_20260928.json'), result)
print(json.dumps(result))
