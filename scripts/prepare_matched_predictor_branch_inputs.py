#!/usr/bin/env python3
"""Prepare all 125 marker slots on all 70 views for matched predictor fits."""
import argparse
import fcntl
import gzip
import json
from pathlib import Path
import shutil

from matched_predictor_branch_inputs import (SCHEMA, load, joint_rows, case, topology,
                                           fasta_text, summary, verify)
from audit_selected_taxon_identity_snapshot_v2 import sha

STATUS = 'complete_full_matched_predictor_branch_inputs_pending_independent_readback'


def run(path):
    path = Path(path); plan = json.loads(path.read_text()); verify(plan['pins'])
    assert plan['schema'] == SCHEMA
    assert shutil.disk_usage('.').free >= plan['resources']['minimum_free_disk_gib']*2**30
    root = Path(plan['output']); root.mkdir(parents=True, exist_ok=False)
    lock = (root/'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    source, bindings = load(); bindings.update(plan['pins']); bindings[str(path)] = sha(path)
    assert (len(source['views']),len(source['markers']),len(source['positions'])) == (70,125,526)
    result, artifacts = produce(source, root)
    assert result['comparison_cases'] == 8750 and result['original_internal_branch_cells'] == 4523750
    verify(bindings)
    receipt = dict(status=STATUS,plan_sha256=sha(path),**result,source_hashes=bindings,artifacts=artifacts,
        scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return receipt


def produce(source, root, verbose=True):
    (root/'input_axes.json').write_text(json.dumps(dict(schema=SCHEMA, views=source['views'],
        positions=source['positions'],markers=source['markers'],taxa=source['taxa']),indent=2)+'\n')
    records = []; inputs = {}; cells = []
    (root/'inputs').mkdir()
    with gzip.open(root/'comparison_cases.jsonl.gz','wt') as f:
        for marker in source['markers']:
            pair = source['alignments'][marker]
            if pair is None: columns, data, account = [], {}, []
            else:
                ctx = {t:r for (m,t),r in source['contexts'].items() if m == marker}
                columns,data,account = joint_rows(*pair,ctx,source['required'][marker])
            cells.extend(dict(marker=marker,**r) for r in account)
            for vi in range(len(source['views'])):
                record, selected = case(source, vi, marker, columns, data)
                records.append(record); f.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
                key = record['input_id']
                if key is not None and key not in inputs:
                    folder = root/'inputs'/key; folder.mkdir()
                    for label, values in selected.items(): (folder/(label+'.faa')).write_text(fasta_text(values))
                    tree = topology(record['taxa'],record['internal_branch_mapping'],source['positions'])
                    (folder/'topology.nwk').write_text(tree)
                    config = {k:record[k] for k in ['input_id','marker','taxa','columns','required_observed','alignment_sha256']}
                    config['internal_splits'] = [r['split_mask_hex'] for r in record['internal_branch_mapping']]
                    config['topology_sha256'] = sha(folder/'topology.nwk')
                    config['branch_initialization'] = 'all_edges_0.1_expected_state_substitutions_per_site_start_only'
                    (folder/'config.json').write_text(json.dumps(config,indent=2)+'\n')
                    inputs[key] = record
            if verbose: print('matched_predictor_marker_prepared',marker,'cases',len(records),'unique_inputs',len(inputs),flush=True)
    (root/'joint_taxon_observation_accounting.json').write_text(json.dumps(cells,indent=2)+'\n')
    (root/'input_manifest.json').write_text(json.dumps([dict(input_id=k,marker=r['marker'],taxa=len(r['taxa']),
        columns=len(r['columns']),internal_splits=len(r['internal_branch_mapping'])) for k,r in sorted(inputs.items())],indent=2)+'\n')
    result = summary(source,records,inputs)
    artifacts = {str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file() and p.name != 'stage.lock'}
    return result, artifacts


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
