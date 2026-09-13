import io
import sys
import unittest
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from convert_esmfold_snapshot import cif_text

class CifConversionTests(unittest.TestCase):
    def pdb(self):
        return ''.join(f'ATOM  {i:5d}  {atom:<3s} ALA A{res:4d}    {x:8.3f}{-2.:8.3f}{3.:8.3f}{1.:6.2f}{87.12:6.2f}           C\n'
                       for i,atom,res,x in [(1,'CA',1,1.123),(2,'CB',1,2.345),(3,'CA',2,4.567),(4,'CB',2,5.678)])+'END\n'
    def test_exact_atom_fields_and_sequence_preserved(self):
        d=MMCIF2Dict(io.StringIO(cif_text(self.pdb(),'AA','ESM-test')))
        self.assertEqual(d['_entity_poly.pdbx_seq_one_letter_code_can'],['AA'])
        self.assertEqual(d['_atom_site.Cartn_x'],['1.123','2.345','4.567','5.678'])
        self.assertEqual(d['_atom_site.label_seq_id'],['1','1','2','2'])
        self.assertEqual(d['_atom_site.B_iso_or_equiv'],['87.12']*4)
        self.assertEqual(d['_atom_site.label_atom_id'],['CA','CB','CA','CB'])
    def test_sequence_and_numbering_errors_rejected(self):
        with self.assertRaises(ValueError):cif_text(self.pdb(),'AC','ESM-test')
        with self.assertRaises(ValueError):cif_text(self.pdb().replace('ALA A   2','ALA A   3'),'AA','ESM-test')
    def test_alternate_atom_requires_explicit_handling(self):
        text=self.pdb();text=text[:16]+'A'+text[17:]
        with self.assertRaises(ValueError):cif_text(text,'AA','ESM-test')

if __name__=='__main__':unittest.main()
