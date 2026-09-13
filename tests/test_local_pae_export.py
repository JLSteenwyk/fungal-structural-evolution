import gzip
import hashlib
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from export_local_marker_pae import encode_pae
from retrieve_marker_pae import validate_pae

class LocalPaeTests(unittest.TestCase):
    def model(self):return {'length':2,'sequence_sha256':hashlib.sha256(b'AC').hexdigest()}
    def test_direction_and_precision_preserved(self):
        pae=np.array([[0.,1.234567],[9.876543,0.]],dtype=np.float32)
        raw,compressed=encode_pae('AC',pae,31.,self.model())
        self.assertEqual(gzip.decompress(compressed),raw)
        np.testing.assert_array_equal(validate_pae(raw,2),pae)
        self.assertEqual(encode_pae('AC',pae,31.,self.model())[1],compressed)
    def test_wrong_sequence_rejected(self):
        with self.assertRaises(ValueError):encode_pae('CA',np.zeros((2,2)),31.,self.model())
    def test_shape_nonfinite_and_maximum_errors_rejected(self):
        for matrix,maximum in [(np.zeros((3,3)),31.),(np.full((2,2),np.nan),31.),(np.full((2,2),32.),31.),(np.zeros((2,2)),float('nan')),(np.zeros((2,2)),0.)]:
            with self.assertRaises(ValueError):encode_pae('AC',matrix,maximum,self.model())

if __name__=='__main__':unittest.main()
