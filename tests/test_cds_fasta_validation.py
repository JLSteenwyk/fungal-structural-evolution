import gzip
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from download_coding_sequences import validate_fasta

class CDSValidation(unittest.TestCase):
    def check(self,text):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'cds.fna.gz'
            with gzip.open(path,'wt') as handle:handle.write(text)
            return validate_fasta(path)

    def test_partial_and_unlinked_records_are_retained_as_flags(self):
        row=self.check('>a [protein_id=P1]\nATG\n>b [pseudo=true]\nATNN\n')
        self.assertEqual(row,{'cds_records':2,'nucleotide_bases':7,'records_without_protein_id':1,'non_triplet_length_records':1,'records_with_ambiguous_bases':1})

    def test_duplicate_records_are_not_silently_merged(self):
        with self.assertRaises(ValueError):self.check('>a\nATG\n>a\nATG\n')

    def test_non_dna_content_is_rejected(self):
        with self.assertRaises(ValueError):self.check('>a\nMPEPTIDE\n')

if __name__=='__main__':unittest.main()
