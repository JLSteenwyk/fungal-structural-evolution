import gzip
import hashlib
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from retrieve_matched_models import decode_http_payload, polymer_sequences


class TransportTests(unittest.TestCase):
    def test_gzip_and_identity_produce_identical_content(self):
        content = b'data_model\n_entry.id model\n'
        raw = gzip.compress(content, mtime=0)
        decoded, receipt = decode_http_payload(raw)
        self.assertEqual(decoded, content)
        self.assertEqual(receipt['response_encoding'], 'gzip')
        self.assertEqual(receipt['response_payload_sha256'], hashlib.sha256(raw).hexdigest())
        self.assertEqual(decode_http_payload(content)[0], content)

    def test_literal_sequence_fallback_does_not_guess_modifications(self):
        self.assertEqual(polymer_sequences({'_entity_poly.pdbx_seq_one_letter_code': ['ACD\nEF']}), ['ACDEF'])
        self.assertEqual(polymer_sequences({'_entity_poly.pdbx_seq_one_letter_code': ['A(MSE)C']}), ['A(MSE)C'])
        self.assertEqual(polymer_sequences({'_entity_poly.pdbx_seq_one_letter_code_can': ['ACD'], '_entity_poly.pdbx_seq_one_letter_code': ['EEE']}), ['ACD'])

    def test_corrupt_gzip_rejected(self):
        compressed = gzip.compress(b'data_model\n', mtime=0)
        with self.assertRaises((EOFError, OSError)):
            decode_http_payload(compressed[:-6])


if __name__ == '__main__':
    unittest.main()
