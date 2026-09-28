#!/usr/bin/env python3
import copy
import json
from pathlib import Path
import sys
import unittest
from ancestral_chain_attempt import sha,write_json
from advance_ancestral_categorical_diagnostics import grouped_jobs,ready,process_group
from test_ancestral_state_quartet import StateQuartetTests


class QueueTests(unittest.TestCase):
    def setUp(self):
        self.fixture=StateQuartetTests();self.fixture.setUp();self.root=self.fixture.root

    def tearDown(self):self.fixture.tearDown()

    def test_four_unique_seeds_required(self):
        self.assertEqual(len(grouped_jobs(self.fixture.jobs)),1)
        bad=copy.deepcopy(self.fixture.jobs);bad[-1]['chain']['seed']=bad[0]['chain']['seed']
        with self.assertRaises(AssertionError):grouped_jobs(bad)
        with self.assertRaises(AssertionError):grouped_jobs(bad[:-1])

    def test_incomplete_or_failed_extraction_waits(self):
        self.assertTrue(ready(self.fixture.jobs,self.fixture.state_root))
        path=self.fixture.state_root/'synthetic-chain-3/disposition.json'
        write_json(path,dict(status='state_extraction_attempt_failed_requires_review'))
        self.assertFalse(ready(self.fixture.jobs,self.fixture.state_root))
        path.unlink();self.assertFalse(ready(self.fixture.jobs,self.fixture.state_root))

    def test_worker_provenance_to_both_reports(self):
        mapping=self.root/'mapping.tsv';mapping.write_text('synthetic mapping\n')
        # Rebind the synthetic extraction fixtures after replacing their mapping hash.
        for job in self.fixture.jobs:
            chain=job['chain']['chain_id'];ap=self.root/chain/'attempt-0001-sample-audit.json'
            a=json.loads(ap.read_text());a['mapping_sha256']=sha(mapping);write_json(ap,a)
            dp=self.fixture.state_root/chain/'disposition.json';d=json.loads(dp.read_text())
            rp=Path(d['receipt']);r=json.loads(rp.read_text());r['summaries'][0]['source_audit_sha256']=sha(ap);write_json(rp,r)
            ar=Path(d['attempt_receipt']);r=json.loads(ar.read_text());r['artifacts'][str(rp.relative_to(ar.parent))]=sha(rp);write_json(ar,r)
            d['receipt_sha256']=sha(rp);d['attempt_receipt_sha256']=sha(ar);write_json(dp,d)
        jobs=self.root/'jobs.json';write_json(jobs,self.fixture.jobs)
        pp=self.root/'producer.json';write_json(pp,dict(jobs=str(jobs),output=str(self.root),iterations=1000,mapping=str(mapping),pins={}))
        ep=self.root/'extraction.json';write_json(ep,dict(producer_plan=str(pp),output=str(self.fixture.state_root),pins={}))
        (self.fixture.state_root/'run_plan.json').write_bytes(ep.read_bytes())
        plan=dict(output=str(self.root/'queue'),producer_plan=str(pp),extraction_plan=str(ep),working_directory=str(Path.cwd()),python=sys.executable,
            diagnostic_python=str(Path('SOFTWARE/ancestral-diagnostics-20260927/bin/python').absolute()),per_process_address_space_bytes=16*1024**3,
            per_group_timeout_seconds=120,pins={'/usr/bin/env':sha('/usr/bin/env'),'/usr/bin/prlimit':sha('/usr/bin/prlimit')})
        result=process_group('synthetic-input',plan)
        self.assertEqual(result['status'],'categorical_reports_complete_not_posterior_qualification')
        self.assertEqual(sha(result['receipt']),result['receipt_sha256'])
        again=process_group('synthetic-input',plan)
        self.assertEqual(result,again)


if __name__=='__main__':unittest.main()
