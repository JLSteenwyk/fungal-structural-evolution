#!/usr/bin/env python3
"""Fresh full readback root; bounded parallel cohorts and complete SQL links."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import csv
import fcntl
import gzip
import itertools
import multiprocessing
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile
import time

from full_covariance_qualification_sources import load,jsonl,audit_id,SUMMARY_FIELDS,LINK_EXTRA
from full_entity_operator_sources import MODES
from full_expanded_model_design_sources import SETTING_FIELDS,digest
from parallel_covariance_readback_v2 import verify_cohort
from reference_measurement_union_sources import bind,verify
from ancestral_chain_attempt import sha


STATUS='passed_full_uniform_covariance_process_latent_and_sql_readback_v2'

# The source is loaded and hash-verified once before forking. Workers receive
# cohort payloads only; the inherited source is never sent through IPC.
_SOURCE=None
_PLAN=None
_FAILURE_DIRECTORY=None


def execute_cohort(cohort,designs,records):
    assert _SOURCE is not None and _PLAN is not None and _FAILURE_DIRECTORY is not None
    result=verify_cohort(_SOURCE,_PLAN,cohort,designs,records,_FAILURE_DIRECTORY)
    result['worker_pid']=os.getpid()
    result['execution_backend']='fork_process'
    return result


def caps(plan):
    group=next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    cg=Path('/sys/fs/cgroup')/group.lstrip('/')
    values={k:(cg/k).read_text().strip() for k in ['cpu.max','memory.max','memory.swap.max']}
    assert values=={'cpu.max':str(plan['resources']['cpus']*100000)+' 100000',
                    'memory.max':str(plan['resources']['memory_gib']*2**30),'memory.swap.max':'0'}
    assert all(os.environ.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    return values


def run(path):
    global _SOURCE,_PLAN,_FAILURE_DIRECTORY
    path=Path(path);parallel=json.loads(path.read_text());verify(parallel['pins']);actual=caps(parallel)
    sp=Path(parallel['source_plan']);plan=json.loads(sp.read_text());source,bindings=load(plan,sp);original=dict(bindings)
    source_root=Path(plan['output']);rp=source_root/'receipt.json';receipt=json.loads(rp.read_text())
    assert receipt['status']=='complete_full_uniform_covariance_qualification_pending_independent_readback'
    assert receipt['plan_sha256']==sha(sp) and receipt['source_contract']==source['contract']
    assert receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    assert set(receipt['artifacts'])=={'stage_plan.json','design_covariance_audits.jsonl.gz','setting_audit_links.tsv.gz'}
    bind(bindings,rp)
    for name,h in receipt['artifacts'].items():bind(bindings,source_root/name,h)
    verify(bindings)
    assert json.loads((source_root/'stage_plan.json').read_text())==dict(plan_sha256=sha(sp),schema='expanded-uniform-covariance-qualification-v1')
    for p,h in parallel['pins'].items():bind(bindings,p,h)
    bind(bindings,path)
    root=Path(parallel['output']);assert shutil.disk_usage(root.parent).free>=parallel['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=False);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=dict(source_plan_sha256=sha(sp),parallel_plan_sha256=sha(path),source_producer_receipt_sha256=sha(rp),schema='full-process-covariance-readback-v2')
    (root/'stage_plan.json').write_text(json.dumps(state,indent=2)+'\n');folder=root/'cohorts';folder.mkdir()
    failures=root/'failures';failures.mkdir();_SOURCE=source;_PLAN=plan;_FAILURE_DIRECTORY=failures
    designs=jsonl(source['root']/'unique_designs.jsonl');records=jsonl(source_root/'design_covariance_audits.jsonl.gz')
    started=time.monotonic();workers=parallel['resources']['workers'];capacity=parallel['resources']['maximum_pending_cohorts']
    assert type(workers) is int and 1<=workers<=8 and capacity==2*workers
    pending={};completed=0;manifests={};peak_pending=0
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
        def finish_some():
            nonlocal completed
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                cohort=pending.pop(future);value=future.result()
                assert value['cohort_id']==cohort['cohort_id'] and value['scientific_eligibility'] is False
                assert value['execution_backend']=='fork_process' and value['worker_pid']!=os.getpid()
                fp=folder/(cohort['cohort_id']+'.json');assert not fp.exists()
                fp.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')
                manifests[cohort['cohort_id']]=dict(cohort_id=cohort['cohort_id'],path=str(fp.relative_to(root)),sha256=sha(fp))
                completed+=1;print('parallel_covariance_cohorts',completed,'/',len(source['cohorts']),flush=True)
        for cohort in source['cohorts']:
            while len(pending)>=capacity:finish_some()
            ds=[next(designs) for _ in range(30)]
            rr=[next(records) for _ in range(30*len(MODES)*len(plan['trees']))]
            future=pool.submit(execute_cohort,cohort,ds,rr);pending[future]=cohort
            peak_pending=max(peak_pending,len(pending))
        assert next(designs,None) is next(records,None) is None
        while pending:finish_some()
    assert not list(failures.iterdir()),'Rejected arithmetic cannot emit a completed readback'
    manifest=[manifests[c['cohort_id']] for c in source['cohorts']]
    assert completed==len(manifest)==len(source['cohorts']) and len(manifests)==completed
    assert {p.name for p in folder.iterdir()}=={c['cohort_id']+'.json' for c in source['cohorts']}
    (root/'cohort_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Independently stream the original recipes again for every serialized
    # cohort identity/index, then every one of the original setting links.
    counts=Counter();linkcounts=Counter();audits=design_count=reviews=0;maximum_error=0.
    designs=jsonl(source['root']/'unique_designs.jsonl');records=jsonl(source_root/'design_covariance_audits.jsonl.gz')
    with tempfile.TemporaryDirectory(prefix='parallel-readback-index-',dir=root) as temp:
        db=sqlite3.connect(str(Path(temp)/'audits.sqlite'));db.execute('PRAGMA cache_size=-65536')
        db.execute('CREATE TABLE audits(aid TEXT PRIMARY KEY,did TEXT,mode TEXT,tree TEXT,status TEXT,UNIQUE(did,mode,tree))')
        for cohort,part in zip(source['cohorts'],manifest):
            fp=root/part['path'];assert sha(fp)==part['sha256'];value=json.loads(fp.read_text())
            ds=[next(designs) for _ in range(30)];rr=[next(records) for _ in range(30*len(MODES)*len(plan['trees']))]
            assert value['source_design_batch_sha256']==digest(ds) and value['source_audit_batch_sha256']==digest(rr)
            assert value['cohort_id']==cohort['cohort_id'] and value['records']==cohort['records']
            assert value['design_rows']==30 and value['audit_rows']==len(rr) and value['scientific_eligibility'] is False
            expected=[];local=Counter()
            for d in ds:
                assert d['cohort_id']==cohort['cohort_id']
                for mode,tree in itertools.product(MODES,plan['trees']):
                    r=rr[len(expected)];identifier=audit_id(source['contract'],d['design_id'],mode,tree)
                    assert r['audit_id']==identifier
                    expected.append([identifier,d['design_id'],mode,tree,r['disposition']]);local[r['disposition']]+=1
            assert value['audit_index']==expected and value['audit_status_counts']==dict(local)
            db.executemany('INSERT INTO audits VALUES(?,?,?,?,?)',expected);db.commit()
            counts.update(local);audits+=len(rr);design_count+=30
            reviews+=value['conservative_independent_classification_differences']
            maximum_error=max(maximum_error,value['maximum_absolute_projected_gram_error'])
        assert next(designs,None) is next(records,None) is None
        settings=links=0
        with gzip.open(source['root']/'model_settings.tsv.gz','rt') as f,gzip.open(source_root/'setting_audit_links.tsv.gz','rt') as g:
            original_settings=csv.DictReader(f,delimiter='\t');exported=csv.DictReader(g,delimiter='\t')
            assert original_settings.fieldnames==SETTING_FIELDS and exported.fieldnames==SETTING_FIELDS+LINK_EXTRA
            for row in original_settings:
                for mode,tree in itertools.product(MODES,plan['trees']):
                    record=next(exported);assert {k:record[k] for k in SETTING_FIELDS}==row
                    identifier=audit_id(source['contract'],row['design_id'],mode,tree)
                    found=db.execute('SELECT did,mode,tree,status FROM audits WHERE aid=?',(identifier,)).fetchone()
                    assert found is not None and found[:3]==(row['design_id'],mode,tree)
                    status=found[3];combined=row['disposition'] if row['disposition']!='ready_for_working_covariance_fit' else status
                    assert [record[k] for k in LINK_EXTRA]==[mode,tree,identifier,status,combined]
                    links+=1;linkcounts[combined]+=1
                settings+=1
            assert next(exported,None) is None
        assert db.execute('SELECT COUNT(*) FROM audits').fetchone()[0]==audits;db.close()
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=settings,unique_cohorts=len(manifest),unique_designs=design_count,
        audit_rows=audits,setting_audit_links=links,audit_status_counts=dict(counts),link_status_counts=dict(linkcounts),trees=plan['trees'],loading_modes=MODES)
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS)
    assert settings==source['design_completion']['model_setting_rows'] and design_count==source['design_completion']['unique_designs']
    assert audits==design_count*len(MODES)*len(plan['trees']) and links==settings*len(MODES)*len(plan['trees'])
    artifacts={str(p.relative_to(root)):sha(p) for p in [root/'stage_plan.json',root/'cohort_manifest.json',*sorted(folder.iterdir())]}
    for name,h in artifacts.items():bind(bindings,root/name,h)
    verify(bindings)
    result=dict(status=STATUS,plan_sha256=sha(sp),parallel_plan_sha256=sha(path),producer_receipt_sha256=sha(rp),
        source_contract=source['contract'],**summary,maximum_absolute_projected_gram_error=maximum_error,
        conservative_independent_classification_differences=reviews,source_hashes=bindings,
        parallel_workers=workers,maximum_observed_pending_cohorts=peak_pending,actual_cgroup_limits=actual,
        execution_backend='fork_process',worker_pids=sorted({json.loads((root/p['path']).read_text())['worker_pid'] for p in manifest}),
        failure_capture_enabled=True,original_numerical_arithmetic_unchanged=True,
        scientific_eligibility=False,elapsed_seconds=time.monotonic()-started,scope=parallel['scope'])
    with (root/'readback.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(summary),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args();run(a.plan)
