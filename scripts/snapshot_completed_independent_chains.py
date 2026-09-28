#!/usr/bin/env python3
"""Freeze and verify available terminal chains; do not infer convergence."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import psutil
from ancestral_chain_attempt import sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    files = sorted(args.root.glob('*/attempt-*-sample-audit.json'))
    assert files
    rows = []
    checked = {}
    def verify(path, expected):
        path = Path(path)
        if str(path) not in checked:
            checked[str(path)] = sha(path)
        assert checked[str(path)] == expected, str(path)
    for path in files:
        audit = json.loads(path.read_text())
        assert audit['status'] == 'all_saved_alignments_and_candidate_nodes_checked'
        rp = Path(audit['attempt_receipt']); verify(rp, audit['attempt_receipt_sha256'])
        receipt = json.loads(rp.read_text()); assert receipt['exit_code'] == 0
        cp = path.parent/'configuration.json'; config = json.loads(cp.read_text())
        encoded = json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False)
        assert hashlib.sha256(encoded.encode()).hexdigest() == receipt['configuration_sha256']
        for name,h in config['pins'].items(): verify(name,h)
        for name,h in receipt['artifacts'].items(): verify(rp.parent/name,h)
        verify(audit['scalar_log'],audit['scalar_log_sha256'])
        mapping=Path('results/ancestral/case-local-trees-20260927-v1/ancestral_node_mapping.tsv')
        verify(mapping,audit['mapping_sha256'])
        identity=json.loads((rp.parent/'process.json').read_text())
        try:
            p=psutil.Process(identity['pid'])
            assert p.create_time()!=identity['created'] or p.status()==psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            pass
        samples=audit['candidate_samples']; by_iteration=Counter(s['iteration'] for s in samples)
        assert by_iteration=={i:4 for i in range(0,audit['iterations']+1,10)}
        assert len({s['source_node'] for s in samples})==4
        rows.append(dict(chain=path.parent.name,seed=config['seed'],model=config['model_input_identity'],
            audit=str(path),audit_sha256=sha(path),configuration=str(cp),configuration_sha256=sha(cp),
            iterations=audit['iterations'],saved_alignments=len(by_iteration),candidate_samples=len(samples),
            elapsed_seconds=receipt['elapsed_seconds']))
    groups=Counter(r['model'] for r in rows)
    write_json(args.output,dict(status='verified_frozen_terminal_chain_subset',created_utc=datetime.now(timezone.utc).isoformat(),
        chains=len(rows),saved_alignments=sum(r['saved_alignments'] for r in rows),candidate_samples=sum(r['candidate_samples'] for r in rows),
        model_chain_counts=dict(groups),models_with_four_chains=sum(v==4 for v in groups.values()),
        checked_source_files=len(checked),rows=rows,script_sha256=sha(__file__),
        scope='Only terminal chains with saved-sample integrity audits present at snapshot time. '
        'No chain mixing, stationary posterior, or overall run completion claim.'))
    print(json.dumps(dict(chains=len(rows),model_chain_counts=dict(groups),checked_source_files=len(checked))))


if __name__=='__main__':
    main()
