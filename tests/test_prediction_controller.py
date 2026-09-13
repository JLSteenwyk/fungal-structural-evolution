import os
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from advance_ecology_predictions import completed_chunk,process_identity


class PredictionControllerTests(unittest.TestCase):
    def test_requires_clean_complete_matching_chunk(self):
        valid={'status':'production_chunk_finished','config_sha256':'pinned','interrupted':False,'remaining_eligible':0,'oom_deferred':0}
        self.assertTrue(completed_chunk(valid,'pinned'))
        for field,value in [('status','running'),('config_sha256','other'),('interrupted',True),('remaining_eligible',1),('oom_deferred',1)]:
            self.assertFalse(completed_chunk(valid|{field:value},'pinned'))
        self.assertFalse(completed_chunk({},'pinned'))

    def test_live_process_identity_and_absent_process(self):
        identity=process_identity(os.getpid())
        self.assertIsNotNone(identity)
        self.assertEqual(identity,process_identity(os.getpid()))
        self.assertIsNone(process_identity(-1))


if __name__=='__main__':unittest.main()
