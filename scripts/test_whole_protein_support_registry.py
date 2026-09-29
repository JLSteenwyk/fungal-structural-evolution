#!/usr/bin/env python3
"""Verify that consumers reject pending, changed and mismatched support inputs."""
import json
import tempfile
import unittest
from pathlib import Path
from whole_protein_support_registry import SupportRegistry
from screen_duplication_alignment_reuse import sha


class RegistryGate(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.certificate = self.root/'certificate.json'
        self.certificate.write_text('{}\n')
        self.table = self.root/'resolved_input_support.jsonl'
        self.table.write_text(json.dumps(dict(fit_input_id='fixture-input', input_sha256='fixture-hash',
            geometry_id='fixture-geometry', classification='zero_supported_to_numeric_tolerance',
            original_classification='unresolved_certificate', source_kind='verified_repaired_witness',
            certificate_file=str(self.certificate), certificate_file_sha256=sha(self.certificate)))+'\n')
        counts = dict(inputs=1, geometries=1, settings=2, changed_inputs=1, changed_settings=2)
        receipt = dict(status='complete_resolved_whole_protein_support_pending_readback',
                       artifacts={self.table.name:sha(self.table)}, **counts)
        self.receipt = self.root/'receipt.json'
        self.receipt.write_text(json.dumps(receipt))
        self.proof = self.root/'proof.json'
        self.proof.write_text(json.dumps(dict(status='passed_full_resolved_whole_protein_support_readback',
            source_receipt_sha256=sha(self.receipt), source_hashes={str(self.certificate):sha(self.certificate)}, **counts)))

    def test_checked_input_resolves_with_original_status(self):
        row = SupportRegistry(self.root, self.proof).resolve('fixture-input', 'fixture-hash')
        self.assertTrue(row['zero_reference_supported'])
        self.assertEqual(row['original_classification'], 'unresolved_certificate')

    def test_wrong_input_hash_is_rejected(self):
        registry = SupportRegistry(self.root, self.proof)
        with self.assertRaises(ValueError): registry.resolve('fixture-input', 'different-hash')

    def test_pending_verification_is_rejected(self):
        p = json.loads(self.proof.read_text()); p['status'] = 'pending'
        self.proof.write_text(json.dumps(p))
        with self.assertRaises(ValueError): SupportRegistry(self.root, self.proof)

    def test_changed_certificate_is_rejected(self):
        self.certificate.write_text('{"changed": true}\n')
        with self.assertRaises(ValueError): SupportRegistry(self.root, self.proof)

    def test_changed_table_is_rejected(self):
        self.table.write_text('')
        with self.assertRaises(ValueError): SupportRegistry(self.root, self.proof)

    def test_verification_scope_mismatch_is_rejected(self):
        p = json.loads(self.proof.read_text()); p['inputs'] = 2
        self.proof.write_text(json.dumps(p))
        with self.assertRaises(ValueError): SupportRegistry(self.root, self.proof)


if __name__ == '__main__': unittest.main()
