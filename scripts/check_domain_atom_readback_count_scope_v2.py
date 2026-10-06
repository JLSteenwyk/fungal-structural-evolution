#!/usr/bin/env python3
"""Check count reporting across every producer shard and reject scope changes."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from review_domain_atom_readback_count_scope_v2 import normalize_counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    producer_path = Path('results/domains/domain-coordinates-20261005-v1/receipt.json')
    producer = json.loads(producer_path.read_text())
    expected = producer['counts']
    observed = {k: v for k, v in expected.items() if k != 'rejected'}
    assert expected['rejected'] == 0
    assert normalize_counts(observed, expected) == expected
    assert normalize_counts(expected, expected) == expected
    for proof in producer['proofs']:
        counts = proof['receipt']['counts']
        assert normalize_counts({k: v for k, v in counts.items() if k != 'rejected'}, counts) == counts
    cases = [
        ('missing_models', {k: v for k, v in observed.items() if k != 'models'}, expected),
        ('missing_backbone_counter', {k: v for k, v in observed.items() if k != 'missing_backbone_intervals'}, expected),
        ('changed_intervals', dict(observed, intervals=observed['intervals'] - 1), expected),
        ('changed_exports', dict(observed, exported=observed['exported'] - 1), expected),
        ('foreign_counter', dict(observed, foreign=0), expected),
        ('boolean_counter', dict(observed, missing_backbone_intervals=False), expected),
        ('negative_counter', dict(observed, rejected=-1), expected),
        ('float_counter', dict(observed, rejected=0.0), expected),
        ('missing_nonzero_rejections', observed, dict(expected, rejected=1)),
        ('changed_explicit_rejections', dict(observed, rejected=1), expected),
        ('invalid_producer_counter', observed, dict(expected, rejected=False)),
    ]
    for label, actual, source in cases:
        try:
            normalize_counts(actual, source)
        except ValueError:
            pass
        else:
            raise AssertionError('Accepted invalid counter: ' + label)
    paths = [producer_path, Path(__file__), Path('scripts/review_domain_atom_readback_count_scope_v2.py')]
    result = dict(status='passed_full_producer_count_scope_and_eleven_rejection_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), producer_shards_checked=len(producer['proofs']),
        global_counts=expected, allowed_missing_counter='rejected', allowed_missing_value=0,
        rejection_controls=[case[0] for case in cases],
        source_hashes={str(path): sha(path) for path in paths},
        full_original_atom_reader_complete=False, native_atom_readback_repeated=False,
        scientific_eligibility=False,
        scope='Reporting software checks across all977actual producer shard count records '
              'and global counts. This is not independent atom-readback completion or native '
              'execution closure. Original full reader remains live and unchanged.')
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
