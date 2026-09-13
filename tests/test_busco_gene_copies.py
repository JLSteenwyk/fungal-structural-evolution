import json
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_busco_gene_copies import classify


def decision(gene, kept=True, status='unique_gene'):
    return {'gene_ids_json': json.dumps([gene]), 'status': status,
            'decision': 'longest_per_gene_lexical_tiebreak' if kept else 'alternative_product_retained_in_source'}


class CopyAnnotationTests(unittest.TestCase):
    def test_isoforms_and_separate_loci_remain_distinct(self):
        hits = [['m', 'Duplicated', 'a'], ['m', 'Duplicated', 'b']]
        isoforms = classify(hits, {'a': decision('g'), 'b': decision('g', False)})
        self.assertEqual(isoforms['duplicate_annotation_category'], 'one_annotated_gene_multiple_products')
        self.assertEqual(isoforms['hit_proteins_retained_as_representatives'], 1)
        loci = classify(hits, {'a': decision('g1'), 'b': decision('g2')})
        self.assertEqual(loci['duplicate_annotation_category'], 'multiple_annotated_genes')
        self.assertEqual(loci['hit_proteins_retained_as_representatives'], 2)

    def test_provisional_orf_never_becomes_confirmed_gene_copy(self):
        result = classify([['m', 'Duplicated', 'a'], ['m', 'Duplicated', 'b']],
                          {'a': decision('g'), 'b': dict(decision('orf', status='provisional_orf'), decision='retained_unresolved_gene')})
        self.assertEqual(result['duplicate_annotation_category'], 'unresolved_gene_mapping')
        self.assertEqual(result['distinct_resolved_gene_ids'], 1)

    def test_original_hit_can_be_excluded_without_recalling_missing_busco(self):
        result = classify([['m', 'Complete', 'a']], {'a': decision('g', False)})
        self.assertEqual(result['busco_status'], 'Complete')
        self.assertEqual(result['hit_proteins_retained_as_representatives'], 0)

    def test_repeated_hit_is_rejected(self):
        with self.assertRaises(ValueError):
            classify([['m', 'Duplicated', 'a'], ['m', 'Duplicated', 'a']], {'a': decision('g')})


if __name__ == '__main__':
    unittest.main()
