import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from extract_creolimax_cds import assemble


def part(start,end,phase,strand='+'):
    return {'sequence':'s','start':start,'end':end,'strand':strand,'phase':str(phase)}

class SplicedCDS(unittest.TestCase):
    def test_split_codon_bases_are_retained(self):
        sequence,_=assemble([part(1,4,0),part(8,12,2)],{'s':'ATGANNNAATTT'})
        self.assertEqual(sequence,'ATGAAATTT')

    def test_negative_strand_uses_transcription_order(self):
        sequence,ordered=assemble([part(1,5,2,'-'),part(9,12,0,'-')],{'s':'AAATTNNNTCAT'})
        self.assertEqual(sequence,'ATGAAATTT')
        self.assertEqual(ordered[0]['start'],9)

    def test_inconsistent_phase_is_not_guessed(self):
        with self.assertRaisesRegex(ValueError,'inconsistent_internal_phase'):
            assemble([part(1,4,0),part(8,12,0)],{'s':'ATGANNNAATTT'})

    def test_out_of_bounds_is_not_truncated(self):
        with self.assertRaisesRegex(ValueError,'coordinates_out_of_bounds'):
            assemble([part(1,6,0)],{'s':'ATG'})

if __name__=='__main__':unittest.main()
