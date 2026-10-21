#!/usr/bin/env python3
"""Close a PAE retry by independently reconciling failure identities and matrices."""
import argparse, gzip, hashlib, json, math, sqlite3, tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def key(model): return f"{model['model_id']}-v{model['version']}"

def rows(path):
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if not line.endswith('\n'): raise ValueError('Partial disposition row')
            yield json.loads(line)

def validate_matrix(raw, model):
    payload=json.loads(raw)
    if not isinstance(payload,list) or len(payload)!=1 or not isinstance(payload[0],dict): raise ValueError('Invalid PAE JSON container')
    matrix=payload[0].get('predicted_aligned_error'); maximum=payload[0].get('max_predicted_aligned_error')
    if not isinstance(matrix,list) or len(matrix)!=model['length'] or not isinstance(maximum,(int,float)) or not math.isfinite(maximum) or maximum<=0: raise ValueError('Invalid PAE dimensions or declared maximum')
    n=0
    for row in matrix:
        if not isinstance(row,list) or len(row)!=model['length']: raise ValueError('Non-square PAE matrix')
        for value in row:
            if not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or value>maximum+.51: raise ValueError('Invalid PAE value')
            n+=1
    return n

def verify_success(row):
    if set(row)!={'model','status','receipt_path','receipt_sha256','receipt'} or row['status']!='retrieval_verified': raise ValueError('Invalid successful retry disposition')
    model,receipt=row['model'],row['receipt']
    if receipt.get('status')!='verified' or any(receipt.get(f)!=model.get(f) for f in ('model_id','version','sequence_sha256','length')) or receipt.get('url')!=model.get('pae_url'): raise ValueError('Retry receipt identity differs')
    rp=Path(row['receipt_path']); matrix=Path(receipt['path'])
    if digest(rp)!=row['receipt_sha256'] or json.loads(rp.read_text())!=receipt: raise ValueError('Retry receipt differs from disposition')
    compressed=matrix.read_bytes()
    if len(compressed)!=receipt['compressed_bytes'] or digest(matrix)!=receipt['gzip_sha256']: raise ValueError('Compressed retry PAE differs')
    raw=gzip.decompress(compressed)
    if len(raw)!=receipt['json_bytes'] or hashlib.sha256(raw).hexdigest()!=receipt['json_sha256']: raise ValueError('Retry JSON differs')
    return receipt['compressed_bytes'],receipt['json_bytes'],validate_matrix(raw,model)

def audit(queue, preserved, recovery, retry_manifest, recovery_readback, retry_receipt, output, receipt, script):
    if output.exists() or receipt.exists(): raise FileExistsError('Fresh output and receipt paths required')
    base=json.loads(recovery_readback.read_text()); retry=json.loads(retry_receipt.read_text())
    if base.get('status')!='passed_order_independent_full_atlas_pae_recovery_readback': raise ValueError('Completed independent recovery readback required')
    if retry.get('status')!='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback': raise ValueError('Completed retry producer receipt required')
    if digest(retry_manifest)!=retry['retry_manifest_sha256']: raise ValueError('Retry manifest digest differs')
    output.mkdir(parents=True); db=sqlite3.connect(output/'retry-source-index.sqlite3')
    try:
        db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
        db.execute('CREATE TABLE source(k TEXT PRIMARY KEY, model TEXT NOT NULL, state TEXT NOT NULL)')
        total=0
        with queue.open() as f:
            for raw in f:
                model=json.loads(raw); db.execute('INSERT INTO source VALUES(?,?,?)',(key(model),json.dumps(model,separators=(',',':')),'unseen')); total+=1
                if total%10000==0: db.commit()
        db.commit(); retained_failures=Counter()
        for manifest,label in ((preserved,'preserved'),(recovery,'recovery')):
            for row in rows(manifest):
                model,status=row.get('model'),row.get('status')
                found=db.execute('SELECT model,state FROM source WHERE k=?',(key(model),)).fetchone()
                if not isinstance(model,dict) or found is None or json.loads(found[0])!=model or found[1]!='unseen': raise ValueError('Base disposition differs from source queue')
                if status=='retrieval_failed': retained_failures[key(model)]=model; state='failed'
                elif status=='retrieval_verified': state='verified'
                else: raise ValueError('Invalid base disposition status')
                db.execute('UPDATE source SET state=? WHERE k=?',(state,key(model)))
        db.commit()
        if db.execute("SELECT COUNT(*) FROM source WHERE state='unseen'").fetchone()[0]: raise ValueError('Base manifests do not cover source queue')
        retry_rows=0; successes=Counter(); errors=Counter(); compressed=json_bytes=cells=0; seen=set()
        for row in rows(retry_manifest):
            model,status=row.get('model'),row.get('status'); k=key(model)
            if k in seen or k not in retained_failures or retained_failures[k]!=model: raise ValueError('Retry disposition is not an exact retained base failure')
            seen.add(k); retry_rows+=1
            if status=='retrieval_verified':
                c,j,n=verify_success(row); compressed+=c;json_bytes+=j;cells+=n;successes['verified']+=1
            elif status=='retrieval_failed' and set(row)=={'model','status','error_type','error'} and isinstance(row['error_type'],str):
                successes['failed']+=1;errors[row['error_type']]+=1
            else: raise ValueError('Invalid retry disposition schema')
        if set(retained_failures)!=seen: raise ValueError('Retry manifest omits or adds retained failure identities')
        if retry_rows != retry['expected_models'] or successes['verified'] != retry['counts']['verified'] or successes['failed'] != retry['counts']['failed']: raise ValueError('Retry producer census differs')
        initial_failed=sum(v for k,v in base['dispositions'].items() if k.endswith('_retrieval_failed'))
        if initial_failed!=retry_rows: raise ValueError('Independent base failed count differs from retry identity census')
        result=dict(status='passed_combined_full_atlas_pae_retry_closure',checked_utc=datetime.now(timezone.utc).isoformat(),expected_models=total,
                    retained_base_failures=retry_rows,retry_dispositions=dict(successes),retry_error_types=dict(errors),
                    final_explicit_failures=successes['failed'],retry_verified_compressed_bytes=compressed,retry_verified_json_bytes=json_bytes,retry_verified_directional_entries=cells,
                    source_queue_sha256=digest(queue),preserved_manifest_sha256=digest(preserved),recovery_manifest_sha256=digest(recovery),retry_manifest_sha256=digest(retry_manifest),
                    recovery_readback_sha256=digest(recovery_readback),retry_receipt_sha256=digest(retry_receipt),independent_reader_sha256=digest(script),scientific_eligibility=False,gpu=False,new_predictions=0,
                    scope='Independent closure: completed base reader is required; this reader independently proves the retry queue is exactly the original retained-failure identity set and validates every successful retry matrix. Explicit remaining failures are retained. This is not confidence calibration, clustering, homology or evolutionary inference.')
        receipt.write_text(json.dumps(result,indent=2)+'\n'); (output/'state.json').write_text(json.dumps(dict(stage='complete',retry_dispositions=dict(successes)),indent=2)+'\n'); return result
    finally: db.close()

def self_test(script):
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary); queue=root/'queue.jsonl'; old=root/'old.jsonl.gz'; recovery=root/'recovery.jsonl.gz'; retry_manifest=root/'retry.jsonl.gz'
        model=dict(model_id='M1',version=6,sequence_sha256='abc',length=2,pae_url='https://x/M1')
        queue.write_text(json.dumps(model)+'\n')
        failed=dict(model=model,status='retrieval_failed',error_type='HTTPError',error='429')
        for path,records in ((old,[failed]),(recovery,[]),(retry_manifest,[failed])):
            with gzip.open(path,'wt') as h:
                for row in records:h.write(json.dumps(row)+'\n')
        base=root/'base.json';base.write_text(json.dumps(dict(status='passed_order_independent_full_atlas_pae_recovery_readback',dispositions=dict(preserved_retrieval_failed=1,recovery_retrieval_failed=0))))
        producer=root/'producer.json';producer.write_text(json.dumps(dict(status='completed_full_atlas_pae_explicit_failure_retry_pending_independent_readback',retry_manifest_sha256=digest(retry_manifest),expected_models=1,counts=dict(attempted=1,verified=0,failed=1))))
        result=audit(queue,old,recovery,retry_manifest,base,producer,root/'out',root/'receipt.json',script)
        assert result['final_explicit_failures']==1 and result['retained_base_failures']==1
    return 'complete-source-failure-set-reconciliation'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('queue','preserved-manifest','recovery-manifest','retry-manifest','recovery-readback','retry-receipt','output','receipt'): p.add_argument('--'+name,type=Path)
    p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test: print(json.dumps(dict(status='passed_pae_retry_reader_literal_control',check=self_test(Path(__file__))),indent=2));return
    if not all([a.queue,a.preserved_manifest,a.recovery_manifest,a.retry_manifest,a.recovery_readback,a.retry_receipt,a.output,a.receipt]): p.error('all paths are required')
    print(json.dumps(audit(a.queue,a.preserved_manifest,a.recovery_manifest,a.retry_manifest,a.recovery_readback,a.retry_receipt,a.output,a.receipt,Path(__file__)),indent=2))
if __name__=='__main__': main()
