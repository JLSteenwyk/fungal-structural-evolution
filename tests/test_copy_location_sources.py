import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_identical_copy_locations import locations

class AnnotationCoordinates(unittest.TestCase):
    def test_creolimax_cds_ids_override_inconsistent_gene_feature_label(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'source.gtf'
            path.write_text('scaffold1\tAUGUSTUS\tgene\t10\t90\t.\t-\t.\tCFRCFRG1\n'
                            'scaffold1\tAUGUSTUS\tCDS\t12\t30\t.\t-\t0\tgene_id "CFRG1"; transcript_id "t1";\n'
                            'scaffold1\tAUGUSTUS\tCDS\t61\t87\t.\t-\t0\tgene_id "CFRG1"; transcript_id "t1";\n')
            self.assertEqual(locations(path, 'creolimax_gtf', {'CFRG1'}), {'CFRG1': {('scaffold1', 12, 87, '-')}})

    def test_multisequence_cds_gene_is_not_merged_into_fake_interval(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'source.gtf'
            path.write_text('s1\ta\tCDS\t12\t30\t.\t+\t0\tgene_id "g";\n'
                            's2\ta\tCDS\t61\t87\t.\t+\t0\tgene_id "g";\n')
            self.assertEqual(len(locations(path, 'creolimax_gtf', {'g'})['g']), 2)

if __name__ == '__main__': unittest.main()
