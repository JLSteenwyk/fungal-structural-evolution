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


class MarkerIntegrity(unittest.TestCase):
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
