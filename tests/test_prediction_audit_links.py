import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_local_predictions import normalize_link

class LinkIdentityTests(unittest.TestCase):
    def test_both_source_link_formats(self):
        digest='a'*64
        for row in ({'sequence_id':'S'+digest},{'sequence_sha256':digest},{'sequence_id':'S'+digest,'sequence_sha256':digest}):
            self.assertEqual(normalize_link(row)['sequence_id'],'S'+digest)
    def test_conflicting_id_rejected(self):
        with self.assertRaises(ValueError):
            normalize_link({'sequence_id':'S'+'a'*64,'sequence_sha256':'b'*64})
    def test_missing_and_malformed_id_rejected(self):
        for row in ({},{'sequence_id':'Sbad'},{'sequence_sha256':'z'*64}):
            with self.assertRaises(ValueError):normalize_link(row)

if __name__=='__main__':unittest.main()
