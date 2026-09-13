import contextlib
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_local_predictions import main


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path,row):path.write_text(json.dumps(row))

class PartialSnapshotTests(unittest.TestCase):
    def fixture(self,root):
        inputs=root/'inputs';pred=root/'pred';inputs.mkdir();pred.mkdir()
        sid='S'+hashlib.sha256(b'AAAA').hexdigest()
        (inputs/'candidates.faa').write_text('>'+sid+'\nAAAA\n')
        (inputs/'all_marker_links.tsv').write_text('sequence_id\tmarker\ttaxon_id\n'+sid+'\tmarker1\tF1\n')
        write_json(inputs/'receipt.json',{'artifacts':{p.name:digest(p) for p in inputs.iterdir()}})
        write_json(pred/'config.json',{'input_receipt_sha256':digest(inputs/'receipt.json'),'max_length':512})
        cp=digest(pred/'config.json')
        # A historical chunk receipt exists, but does not describe the current file set.
        write_json(pred/'last_chunk.json',{'config_sha256':cp,'cached_predictions':0,'new_predictions':128,'remaining_eligible':1000,'length_or_alphabet_deferred':0})
        pdb=pred/(sid+'.pdb')
        pdb.write_text(''.join(f'ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{80.:6.2f}           C\n' for i in range(1,5))+'END\n')
        npz=pred/(sid+'.npz');np.savez(npz,sequence=np.array('AAAA'),ca_plddt=np.full(4,80.),pae=np.zeros((4,4)),max_pae=np.array(31.))
        write_json(pred/(sid+'.json'),{'status':'verified_prediction','sequence_id':sid,'sequence_sha256':sid[1:],'max_predicted_aligned_error':31.,'config_sha256':cp,'length':4,'mean_ca_plddt':80.,'fraction_ca_plddt_below50':0.,'inference_seconds':1.,'peak_gpu_allocated_bytes':10,'artifacts':{p.name:digest(p) for p in [pdb,npz]}})
        return inputs,pred

    def test_partial_snapshot_does_not_claim_historical_chunk_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);inputs,pred=self.fixture(root);out=root/'audit'
            argv=['audit','--inputs',str(inputs),'--predictions',str(pred),'--output',str(out),'--snapshot-live']
            with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):main()
            result=json.loads((out/'receipt.json').read_text())
            self.assertEqual(result['status'],'complete_artifact_readback_of_partial_prediction_snapshot')
            self.assertEqual(result['predictions'],1)
            self.assertIsNone(result['chunk_receipt_sha256'])
            self.assertIsNone(result['remaining_eligible'])

    def test_completed_chunk_mode_still_rejects_count_disagreement(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);inputs,pred=self.fixture(root)
            argv=['audit','--inputs',str(inputs),'--predictions',str(pred),'--output',str(root/'audit')]
            with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'Completed chunk'):main()

    def test_inconsistent_prediction_metadata_rejected(self):
        for field, value in [('sequence_sha256', 'a'*64), ('length', 5),
                             ('fraction_ca_plddt_below50', .25),
                             ('max_predicted_aligned_error', 32.)]:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                root=Path(directory);inputs,pred=self.fixture(root)
                path=next(pred.glob('S*.json'));row=json.loads(path.read_text())
                row[field]=value;write_json(path,row)
                argv=['audit','--inputs',str(inputs),'--predictions',str(pred),
                      '--output',str(root/'audit'),'--snapshot-live']
                with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(ValueError):main()

    def test_complete_grid_and_wrong_remaining_count(self):
        for remaining in [0, 1]:
            with self.subTest(remaining=remaining), tempfile.TemporaryDirectory() as directory:
                root=Path(directory);inputs,pred=self.fixture(root)
                write_json(pred/'last_chunk.json',{
                    'status':'production_chunk_finished',
                    'config_sha256':digest(pred/'config.json'),
                    'cached_predictions':0,'new_predictions':1,
                    'remaining_eligible':remaining,'length_or_alphabet_deferred':0})
                argv=['audit','--inputs',str(inputs),'--predictions',str(pred),
                      '--output',str(root/'audit')]
                with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):
                    if remaining:
                        with self.assertRaisesRegex(ValueError,'identity grid'):main()
                    else:
                        main()
                        self.assertEqual(json.loads((root/'audit/receipt.json').read_text())['remaining_eligible'],0)

if __name__=='__main__':unittest.main()
