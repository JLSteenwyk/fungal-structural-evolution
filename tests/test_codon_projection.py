import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from project_published_outgroup_codons import project

class AnnotatedCodonSpan(unittest.TestCase):
    def test_annotated_initial_and_terminal_partial_bases(self):
        parts=[{'sequence':'s','start':1,'end':5,'strand':'+','phase':'2'},
               {'sequence':'s','start':9,'end':12,'strand':'+','phase':'0'}]
        dna,initial,tail,_=project(parts,{'s':'CCATGNNNAAAT'},'CCATGAAAT')
        self.assertEqual((dna,initial,tail),('ATGAAA',2,1))

    def test_reverse_strand_partial_boundary(self):
        parts=[{'sequence':'s','start':1,'end':9,'strand':'-','phase':'2'}]
        dna,initial,tail,_=project(parts,{'s':'ATTTCATGG'},'CCATGAAAT')
        self.assertEqual((dna,initial,tail),('ATGAAA',2,1))

    def test_genome_disagreement_is_not_overridden_by_protein_match(self):
        parts=[{'sequence':'s','start':1,'end':6,'strand':'+','phase':'0'}]
        with self.assertRaisesRegex(ValueError,'genome_annotation_differs'):
            project(parts,{'s':'ATGAAA'},'ATGAAG')

if __name__=='__main__':unittest.main()
