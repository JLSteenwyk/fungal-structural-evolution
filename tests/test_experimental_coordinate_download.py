import gzip,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from retrieve_experimental_coordinates import verify_gzip


class CoordinateDownloadTests(unittest.TestCase):
    def test_identity_and_complete_stream(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.gz';data=b'data_1ABC\n#\n'+b'x'*2000000;p.write_bytes(gzip.compress(data))
            self.assertEqual(verify_gzip(p,'1ABC'),len(data))
            with self.assertRaises(ValueError):verify_gzip(p,'2ABC')
    def test_truncated_archive_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.gz';p.write_bytes(gzip.compress(b'data_1ABC\n#\n')[:-4])
            with self.assertRaises((EOFError,OSError)):verify_gzip(p,'1ABC')


if __name__=='__main__':unittest.main()
