"""Scientific data-integrity checks for marker extraction and alignment."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


prepare = module('prepare_phylogenetic_markers')
align = module('align_phylogenetic_markers')
assessment = module('assess_marker_alignments')
matrix_builder = module('build_species_matrix')
gene_trees = module('run_marker_gene_trees')


class MarkerIntegrity(unittest.TestCase):
    def test_gene_tree_coverage_filter_counts_only_observed_amino_acids(self):
        kept, removed, threshold = gene_trees.select_rows({
            'A': 'A' * 60 + '-' * 140, 'B': 'A' * 59 + 'X' * 141, 'C': '-' * 200})
        self.assertEqual(threshold, 60)
        self.assertEqual(set(kept), {'A'})
        self.assertEqual(removed, {'B': 59, 'C': 0})

    def test_concatenation_preserves_missingness_and_coordinates(self):
        matrix, partitions, sites, coverage = matrix_builder.concatenate(
            ['A', 'B', 'C'], [('m1', {'A': 'ACX', 'B': '---'}, [1, 3]),
                            ('m2', {'B': 'DF', 'C': 'EF'}, [1, 2])])
        self.assertEqual(matrix, {'A': 'AX--', 'B': '--DF', 'C': '--EF'})
        self.assertEqual(partitions, [('m1', 1, 2), ('m2', 3, 4)])
        self.assertEqual(sites, [(1, 'm1', 1), (2, 'm1', 3), (3, 'm2', 1), (4, 'm2', 2)])
        self.assertEqual(coverage, {'A': 1, 'B': 2, 'C': 2})
        with self.assertRaises(ValueError):
            matrix_builder.concatenate(['A', 'B'], [('m', {'A': 'AC'}, [1, 2])])

    def test_occupancy_masks_and_informative_sites(self):
        statistics, masks = assessment.assess(['AAX-', 'ACX-', 'DC--', 'DC--'])
        self.assertEqual(statistics['parsimony_informative_columns'], 1)
        self.assertEqual(statistics['all_gap_or_ambiguous_columns'], 2)
        self.assertEqual(masks['0.5'], [1, 2])
        self.assertEqual(masks['0.75'], [1, 2])
        with self.assertRaises(ValueError):
            assessment.assess(['AA', 'A'])

    def test_alignment_preserves_residues_and_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / 'input', Path(tmp) / 'output'
            source.write_text('>A\nACDE\n>B\nADE\n')
            target.write_text('>B\nA-DE\n>A\nACDE\n')
            self.assertEqual(align.validate_alignment(source, target), (2, 4))
            for invalid in ['>A\nACDE\n>B\nACDE\n', '>A\nACDE\n>A\nACDE\n', '>A\nACDE\n>B\nADE\n']:
                target.write_text(invalid)
                with self.assertRaises(ValueError):
                    align.validate_alignment(source, target)

    def test_missing_taxon_blocks_complete_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'metadata').mkdir()
            (root / 'metadata/analysis_manifest.tsv').write_text('taxon_id\tstudy_role\nA\tingroup\n')
            (root / 'metadata/qc_input_receipts.json').write_text('[]')
            with self.assertRaisesRegex(ValueError, 'await successful QC'):
                prepare.prepare(root, root / 'out')
            self.assertFalse((root / 'out').exists())

    def test_source_sequence_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'metadata').mkdir()
            (root / 'metadata/analysis_manifest.tsv').write_text('taxon_id\tstudy_role\nA\tingroup\n')
            source = root / 'input.faa'
            source.write_text('>protein\nACDE\n')
            digest = prepare.sha(source)
            (root / 'metadata/qc_input_receipts.json').write_text(json.dumps([
                {'taxon_id': 'A', 'sha256': digest, 'input_path': 'input.faa'}]))
            run = root / 'results/busco/A/run_eukaryota_odb12.2'
            seqdir = run / 'busco_sequences/single_copy_busco_sequences'
            seqdir.mkdir(parents=True)
            (root / 'results/busco/A.receipt.json').write_text(json.dumps({
                'returncode': 0, 'input_sha256': digest}))
            (run / 'full_table.tsv').write_text('m0\tComplete\tprotein\n' + ''.join(
                f'm{i}\tMissing\n' for i in range(1, 125)))
            (seqdir / 'm0.faa').write_text('>protein\nAAAA\n')
            with self.assertRaisesRegex(ValueError, 'differs from source'):
                prepare.prepare(root, root / 'out')
            self.assertFalse((root / 'out').exists())
            (seqdir / 'm0.faa').write_text('>protein\nACDE\n')
            result = prepare.prepare(root, root / 'out')
            self.assertEqual(result['status'], 'complete_extraction')
            self.assertEqual(result['sequences'], 1)
            self.assertEqual((root / 'out/unaligned/m0.faa').read_text(), '>A\nACDE\n')
            source.write_text('>protein\nAAAA\n')
            with self.assertRaisesRegex(ValueError, 'changed BUSCO input'):
                prepare.prepare(root, root / 'other')


if __name__ == '__main__':
    unittest.main()
