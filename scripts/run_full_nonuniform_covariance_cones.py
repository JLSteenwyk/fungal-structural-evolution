#!/usr/bin/env python3
"""Export/read back every exact cone while keeping residual D separate from I."""
import argparse
from collections import Counter
import fcntl
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from full_nonuniform_covariance_cone_sources import load
from independent_nonuniform_covariance_cone import readback
from nonuniform_covariance_cone import record
from reference_measurement_union_sources import bind, verify

PRODUCER = 'complete_full_positive_diagonal_covariance_cones_pending_independent_readback_v1'
READER = 'passed_full_positive_diagonal_covariance_cones_fraction_readback_v1'
SUMMARY = ['logical_cases', 'cohorts', 'certificates', 'case_row_occurrences', 'retained_basis_counts',
           'loading_modes', 'residual_diagonal_prepared', 'raw_reml_basis_qualification_complete',
           'nonuniform_weighting_accepted']


def run(path, reader=False):
    plan = json.loads(path.read_text()); parents, closed, contract, bindings = load(plan, path)
    root = Path(plan['output'])
    assert shutil.disk_usage(root.parent).free >= plan['resources']['minimum_free_disk_gib'] * 2**30
    count = Counter(); occurrences = 0
    summary = dict(logical_cases=closed['logical_cases'], cohorts=closed['cohorts'], certificates=len(parents),
                   loading_modes=['signed', 'unsigned'], residual_diagonal_prepared=False,
                   raw_reml_basis_qualification_complete=False, nonuniform_weighting_accepted=False)
    data = root / 'cohort_cones.jsonl'
    if reader:
        rp = root / 'receipt.json'; producer = json.loads(rp.read_text())
        assert producer['status'] == PRODUCER and producer['plan_sha256'] == sha(path)
        assert producer['source_contract'] == contract and producer['source_hashes'] == bindings
        assert producer['scientific_eligibility'] is False
        assert producer['artifacts'] == {'cohort_cones.jsonl': sha(data)}
        bind(bindings, rp); bind(bindings, data)
        assert not (root / 'readback.json').exists()
        with data.open() as handle:
            for parent in parents:
                saved = json.loads(next(handle)); readback(saved, parent, contract)
                count[str(len(saved['variance_map']['retained_names']))] += 1
                occurrences += saved['records']
            assert next(handle, None) is None
    else:
        root.mkdir(exist_ok=False)
        lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with data.open('x') as handle:
            for parent in parents:
                saved = record(parent, contract)
                handle.write(json.dumps(saved, sort_keys=True, allow_nan=False) + '\n')
                count[str(len(saved['variance_map']['retained_names']))] += 1
                occurrences += saved['records']
    summary.update(case_row_occurrences=occurrences, retained_basis_counts=dict(count))
    if reader:
        assert all(producer[k] == v for k, v in summary.items())
    verify(bindings)
    result = dict(status=READER if reader else PRODUCER, plan_sha256=sha(path), source_contract=contract,
                  **summary, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    if reader:
        result['producer_receipt_sha256'] = sha(rp)
    else:
        result['artifacts'] = {'cohort_cones.jsonl': sha(data)}
    output = root / ('readback.json' if reader else 'receipt.json')
    with output.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(dict(status=result['status'], **summary)), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--reader', action='store_true')
    args = parser.parse_args(); run(args.plan, args.reader)
