import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from link_functional_sites_structures import aggregate_sites


class FunctionalSiteJoinTests(unittest.TestCase):
    def row(self,**kw):
        r={'hit_id':'H1','profile_match_position':'4','sequence_id':'S1','pfam_accession':'PF1.1','pfam_name':'family','pfam_type':'Domain','protein_residue_1based':'10','observed_amino_acid':'A','overlaps_another_ga_hit':'False','hmm_coverage':'1','expected_amino_acid':'C','reference_uniprot_accession':'P1','pattern_status':'pattern_residues_not_all_conserved','candidate_active_site':'False'}
        r.update(kw);return r
    def test_nonconserved_and_gap_rows_retained(self):
        result=aggregate_sites([self.row(),self.row(hit_id='H2',protein_residue_1based='',observed_amino_acid='-',pattern_status='pattern_has_unaligned_residue')])
        self.assertEqual(len(result),2)
        self.assertFalse(any(r['conserved_candidate'] for r in result))
    def test_reference_duplicates_not_independent_sites(self):
        r=aggregate_sites([self.row(),self.row(reference_uniprot_accession='P2')])
        self.assertEqual(len(r),1);self.assertEqual(r[0]['reference_accessions'],'P1;P2')
    def test_inconsistent_projection_rejected(self):
        with self.assertRaises(ValueError):aggregate_sites([self.row(),self.row(protein_residue_1based='11')])


if __name__=='__main__':unittest.main()
