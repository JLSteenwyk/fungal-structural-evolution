#!/usr/bin/env python3
"""Retrieve the complete independently reconstructed missing AFDB PAE queue with two HTTP workers."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import shutil
import time

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify
from retrieve_marker_pae import retrieve


def completed_requests(models, executor, submit, window):
    """Keep a fixed future window rather than materializing millions of futures."""
    models = iter(models); pending = {}; exhausted = False
    while pending or not exhausted:
        while not exhausted and len(pending) < window:
            try:model=next(models)
            except StopIteration:exhausted=True;break
            pending[executor.submit(submit,model)]=model
        if not pending:break
        ready,_=wait(pending,return_when=FIRST_COMPLETED)
        for future in ready:
            model=pending.pop(future)
            yield model,future


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.receipt.exists()
    plan=json.loads(a.plan.read_text());pins=dict(plan['pins']);verify(pins)
    audit_path,tp=Path(plan['availability_readback']),Path(plan['availability_transport'])
    audit,t=[json.loads(path.read_text()) for path in (audit_path,tp)]
    assert audit['status']=='passed_full_atlas_pae_model_and_missing_queue_reconstruction'
    assert audit['missing_retrievable_afdb_models']==plan['expected_models']
    assert t['status']=='verified_original_software_wait_and_whole_wrapper_payloads'
    assert t['validation_sha256']==sha(audit_path) and t['original_tool_terminal_exit_code']==0
    assert t['whole_wrapper_initial_and_terminal_payloads_matched']
    assert t['manager_start_records']==t['manager_completion_records']==1
    verify(pins)
    root,cache=Path(plan['output']),Path(plan['cache']).resolve();assert not root.exists() and not cache.exists()
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(parents=True);cache.mkdir(parents=True)
    counts=Counter(attempted=0,verified=0,failed=0,compressed_bytes=0,json_bytes=0)
    errors=Counter();started=time.monotonic();manifest=root/'retrieval_dispositions.jsonl.gz'
    def fetch(model):
        if shutil.disk_usage(cache).free<plan['resources']['emergency_free_disk_gib']*2**30:
            raise RuntimeError('Emergency disk reserve reached')
        directory=cache/model['sequence_sha256'][:2];directory.mkdir(exist_ok=True)
        # Original qualified URL/provenance/hash/shape/finite/max-bound checks,
        # including its four bounded HTTP attempts and original rounding bound.
        return retrieve(model,directory)
    def state(status):
        value=dict(stage=status,counts=dict(counts),error_types=dict(errors),expected_models=plan['expected_models'],
                   elapsed_seconds=time.monotonic()-started,http_workers=2,future_window=plan['future_window'])
        temporary=root/'state.tmp';temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(root/'state.json')
    state('retrieving_complete_missing_queue')
    raw=manifest.open('xb');gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)
    try:
        with io.TextIOWrapper(gz,encoding='utf-8',newline='') as out,Path(plan['queue']).open() as queue,ThreadPoolExecutor(max_workers=2) as executor:
            models=(json.loads(line) for line in queue)
            for model,future in completed_requests(models,executor,fetch,plan['future_window']):
                counts['attempted']+=1
                try:r=future.result()
                except Exception as error:
                    if isinstance(error,RuntimeError) and str(error)=='Emergency disk reserve reached':raise
                    errors[type(error).__name__]+=1;counts['failed']+=1
                    row=dict(model=model,status='retrieval_failed',error_type=type(error).__name__,error=str(error))
                else:
                    assert r['status']=='verified'
                    for key in ('model_id','version','sequence_sha256','length'):assert r[key]==model[key]
                    rp=Path(r['path']).with_name(f"{model['model_id']}-v{model['version']}.receipt.json")
                    assert json.loads(rp.read_text())==r
                    counts['verified']+=1;counts['compressed_bytes']+=r['compressed_bytes'];counts['json_bytes']+=r['json_bytes']
                    row=dict(model=model,status='retrieval_verified',receipt_path=str(rp),receipt_sha256=sha(rp),receipt=r)
                out.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
                if counts['attempted']%100==0:
                    out.flush();state('retrieving_complete_missing_queue')
                if counts['compressed_bytes']>plan['resources']['maximum_cache_bytes']:
                    raise RuntimeError('Declared cache output allowance reached')
    finally:raw.close()
    assert counts['attempted']==plan['expected_models']==counts['verified']+counts['failed']
    verify(pins)
    for path in (a.plan,audit_path,tp,manifest,Path(__file__)):bind(pins,path)
    # Individual receipt/matrix checksums remain inside the bound complete gzip
    # disposition archive; do not create a millions-entry public JSON pin map.
    verify(pins)
    result=dict(status='completed_complete_missing_afdb_pae_retrieval_attempts_pending_independent_readback',
                checked_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),error_types=dict(errors),
                source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
                scope='Every independently reconstructed missing AFDB model attempted with original '
                      'qualified version/sequence/URL and matrix checks, two HTTP workers and bounded futures. '
                      'All successes and errors retained with source model and per-file checksum provenance. '
                      'Attempt completion is not universal availability; full independent source/receipt/matrix '
                      'readback, fresh union, native context confidence/calibration and biological analyses remain required.')
    with a.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    state('all_attempts_complete_pending_original_wait_and_independent_readback')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
