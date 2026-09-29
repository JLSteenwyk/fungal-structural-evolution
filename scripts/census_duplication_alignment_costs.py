#!/usr/bin/env python3
"""Describe recorded alignment runtimes; these are not prospective ETAs."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt_path = args.source / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'complete_duplication_alignment_dispositions_pending_readback'
    assert sha(args.plan) == receipt['plan_sha256']
    bindings = {str(receipt_path): sha(receipt_path), str(args.plan): sha(args.plan),
                __file__: sha(__file__)}
    manifest = args.source / 'checkpoint_manifest.tsv'
    bindings[str(manifest)] = receipt['artifacts']['checkpoint_manifest.tsv']
    assert sha(manifest) == bindings[str(manifest)]
    groups = defaultdict(list)
    counts = Counter()
    identities = set()
    total_bytes = 0
    args.output.mkdir(parents=True, exist_ok=False)
    with manifest.open() as handle:
        for i, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
            path = args.source / row['path']
            raw = path.read_bytes()
            assert hashlib.sha256(raw).hexdigest() == row['sha256'], path
            record = json.loads(raw)
            identity = (record['pair_key'], record['mask'], record['order'])
            assert identity not in identities
            assert record['mask'] in ('full', 'plddt70') and record['order'] in (0, 1)
            assert record['status'] == row['status']
            assert record['plan_sha256'] == receipt['plan_sha256']
            assert record['input_manifest_sha256'] == receipt['input_manifest_sha256']
            identities.add(identity)
            seconds = record['elapsed_seconds']
            assert math.isfinite(seconds) and seconds >= 0
            if record['status'] == 'aligned':
                length = max(record['metrics']['length_left'], record['metrics']['length_right'])
                assert length > 0
                length_bin = next((str(limit) for limit in (250, 500, 1000, 2000, 4000)
                                   if length <= limit), '>4000')
            else:
                length_bin = 'not_aligned'
            groups[(record['mask'], record['status'], length_bin)].append(seconds)
            counts[record['mask'] + ':' + record['status']] += 1
            total_bytes += len(raw)
            if i % 25000 == 0:
                print('Checked alignment timing records', i, flush=True)
    assert len(identities) == receipt['directed_dispositions']
    assert dict(counts) == receipt['counts']
    pairs = {key[0] for key in identities}
    assert len(pairs) == receipt['distinct_model_pairs']
    assert all((pair, mask, order) in identities for pair in pairs
               for mask in ('full', 'plddt70') for order in (0, 1))
    rows = []
    for (mask, status, length_bin), values in sorted(groups.items()):
        ordered = sorted(values)
        rows.append(dict(mask=mask, status=status, maximum_length_upper_bin=length_bin,
                         records=len(values), summed_elapsed_seconds=math.fsum(values),
                         median_seconds=statistics.median(values),
                         p95_seconds=ordered[math.ceil(.95 * len(values)) - 1],
                         maximum_seconds=max(values)))
    result = dict(status='complete_recorded_duplication_alignment_cost_census',
                  directed_dispositions=len(identities), distinct_pairs=len(pairs),
                  checkpoint_bytes=total_bytes, strata=rows, source_hashes=bindings,
                  scope='Every checkpoint hash, identity, input binding, disposition and recorded elapsed time checked. Length bins use the longer input and inclusive upper limits. Summed elapsed time is concurrent task wall time, not measured CPU time or pipeline wall time. No extrapolation to new pairs, numerical alignment validation, or result reuse authorization.')
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('strata', 'source_hashes')}), flush=True)


if __name__ == '__main__':
    main()
