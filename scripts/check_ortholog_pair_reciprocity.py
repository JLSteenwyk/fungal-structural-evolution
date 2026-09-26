#!/usr/bin/env python3
"""Reproducible independent count checks for the compact ortholog pair stream.

Run: python scripts/check_ortholog_pair_reciprocity.py
Uses temporary files only; never touches the live production outputs.
"""
from collections import Counter
import hashlib
from pathlib import Path
import random
import tempfile
import unittest

from audit_ortholog_pair_reciprocity import encode_pair, summarize_sorted


class PairStreamChecks(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.path = Path(self.scratch.name) / 'pairs.hex'

    def compare(self, pairs, limit=1000):
        # The oracle operates on directed tuples, independently of stream runs.
        directed = Counter(pairs)
        unordered = {tuple(sorted(p)) for p in directed}
        multiplicities = [(a, b, directed[a, b], directed[b, a])
                          for a, b in sorted(unordered)]
        expected = {
            'directed_incidences': len(pairs),
            'unique_unordered_pairs': len(unordered),
            'unique_directed_pairs': len(directed),
            'duplicate_directed_incidences': len(pairs) - len(directed),
            'pairs_with_repeated_direction': sum(f > 1 or r > 1 for a, b, f, r in multiplicities),
            'pairs_missing_reverse': sum(f == 0 or r == 0 for a, b, f, r in multiplicities),
            'pairs_with_unequal_multiplicity': sum(f != r for a, b, f, r in multiplicities),
            'pairs_exactly_once_each_direction': sum(f == r == 1 for a, b, f, r in multiplicities),
            'violation_examples': [dict(left_id=a, right_id=b, forward=f, reverse=r)
                                   for a, b, f, r in multiplicities if (f, r) != (1, 1)][:limit],
            'example_limit': limit,
        }
        raw = b''.join(sorted(encode_pair(a, b) for a, b in pairs))
        self.path.write_bytes(raw)
        expected['sorted_pairs_sha256'] = hashlib.sha256(raw).hexdigest()
        self.assertEqual(summarize_sorted(self.path, len(pairs), limit), expected)

    def test_empty_and_boundary_ids(self):
        self.compare([])
        self.compare([(0, 2**24 - 1), (2**24 - 1, 0)])
        self.assertEqual(encode_pair(0, 2**24 - 1), b'000000ffffff0\n')
        self.assertEqual(encode_pair(2**24 - 1, 0), b'000000ffffff1\n')

    def test_multiplicity_and_sample_limit(self):
        pairs = [(0, 1), (1, 0), (2, 3), (2, 3), (3, 2), (4, 5),
                 (6, 7), (6, 7), (7, 6), (7, 6), (9, 8)]
        for limit in [0, 1, 2, 1000]:
            with self.subTest(limit=limit):
                self.compare(pairs, limit)

    def test_deterministic_random_multisets(self):
        rng = random.Random(20260926)
        for case in range(100):
            pairs = [tuple(rng.sample(range(80), 2)) for _ in range(rng.randrange(300))]
            pairs += pairs[:len(pairs) // 3]
            pairs += [(b, a) for a, b in pairs[:len(pairs) // 2]]
            rng.shuffle(pairs)
            with self.subTest(case=case):
                self.compare(pairs, case % 11)

    def test_reject_invalid_ids(self):
        for pair in [(1, 1), (-1, 2), (0, 2**24), (2**24, 0)]:
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                encode_pair(*pair)

    def test_reject_corrupt_streams(self):
        cases = [b'0000000000010', b'0000000000012\n', b'00000g0000010\n',
                 b'00000000000A0\n', b'0000010000000\n', b'0000010000010\n',
                 b'0000000000010\r\n',
                 b'0000000000020\n0000000000010\n',
                 b'0000000000011\n0000000000010\n']
        for raw in cases:
            with self.subTest(raw=raw):
                self.path.write_bytes(raw)
                with self.assertRaises(ValueError):
                    summarize_sorted(self.path, raw.count(b'\n'))

    def test_reject_count_mismatch(self):
        self.path.write_bytes(encode_pair(0, 1))
        for expected in [0, 2]:
            with self.subTest(expected=expected), self.assertRaises(ValueError):
                summarize_sorted(self.path, expected)


if __name__ == '__main__':
    unittest.main(verbosity=2)
