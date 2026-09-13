import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import match_uniprot_structures as matching

class TaxonomyOverrideTests(unittest.TestCase):
    def test_evidence_checked_and_conflicting_identity_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'config').mkdir()
            source=root/'source.xml'; source.write_text('pinned source')
            evidence=root/'evidence.json'
            evidence.write_text(json.dumps({'status':'verified_deposition_taxonomy_link','taxon_id':'T1','species_taxid':'123','sources':{'s':{'path':'source.xml','sha256':hashlib.sha256(source.read_bytes()).hexdigest()}}}))
            (root/'config/taxonomy_overrides.json').write_text(json.dumps({'T1':{'species_name':'Species one','species_taxid':'123','evidence_path':'evidence.json','evidence_sha256':hashlib.sha256(evidence.read_bytes()).hexdigest()}}))
            with patch.object(matching,'ROOT',root):
                row={'taxon_id':'T1','species_name':'Species one','species_taxid':''}
                self.assertEqual(matching.apply_taxonomy_overrides([dict(row)])[0]['species_taxid'],'123')
                with self.assertRaises(ValueError):matching.apply_taxonomy_overrides([row|{'species_name':'Species two'}])
                source.write_text('changed source')
                with self.assertRaises(ValueError):matching.apply_taxonomy_overrides([dict(row)])

if __name__ == '__main__':
    unittest.main()
