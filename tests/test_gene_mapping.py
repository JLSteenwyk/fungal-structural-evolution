"""Ensure locus mapping does not confuse CDS segments, isoforms and paralogs."""
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('mapping', Path(__file__).resolve().parents[1] / 'scripts/map_proteins_to_genes.py')
mapping = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mapping)


class GeneMapping(unittest.TestCase):
    def test_isoforms_multisegment_and_multiple_loci(self):
        features = [('gene', 'ID=g1'), ('gene', 'ID=g2'),
                    ('mRNA', 'ID=t1;Parent=g1'), ('mRNA', 'ID=t2;Parent=g1'),
                    ('CDS', 'ID=c1;Parent=t1;protein_id=p1'),
                    ('CDS', 'ID=c1;Parent=t1;protein_id=p1'),
                    ('CDS', 'ID=c2;Parent=t2;protein_id=p2;partial=true'),
                    ('CDS', 'ID=c3;Parent=g1,g2;protein_id=p3'),
                    ('CDS', 'ID=c4;Parent=missing;protein_id=p4')]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'x.gff.gz'
            with gzip.open(path, 'wt') as out:
                for kind, attrs in features:
                    out.write(f'ctg\ttest\t{kind}\t1\t10\t.\t+\t.\t{attrs}\n')
            loci, partial = mapping.parse_gff(path)
            self.assertEqual(loci, {'p1': {'g1'}, 'p2': {'g1'}, 'p3': {'g1', 'g2'}, 'p4': set()})
            self.assertEqual(partial, {'p2'})

    def test_escaped_delimiters(self):
        self.assertEqual(mapping.attributes('ID=a%2Cb;Parent=x,y'), {'ID': ['a,b'], 'Parent': ['x', 'y']})

    def test_parent_cycles_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'x.gff.gz'
            with gzip.open(path, 'wt') as out:
                out.write('ctg\ttest\tCDS\t1\t10\t.\t+\t.\tID=c;Parent=c;protein_id=p\n')
            with self.assertRaisesRegex(ValueError, 'Cycle'):
                mapping.parse_gff(path)


if __name__ == '__main__':
    unittest.main()
