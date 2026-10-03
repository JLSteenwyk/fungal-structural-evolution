#!/usr/bin/env python3
"""Checkpoint the complete qualified full-grid fits, retaining every review row."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import sqlite3
import time
from full_retained_shared_entity_fit_sources import load,cohorts,cases,candidate_id,FIT_LINK_EXTRA,SUMMARY_FIELDS,NEW_LINK_EXTRA,SCHEMA,PRODUCER_STATUS,operators_for
from retained_shared_entity_candidate import candidate
from full_expanded_model_design_sources import SETTING_FIELDS
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha


def run(path,stop_after_cohorts=None):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);root=Path(plan['output'])
    assert shutil.disk_usage(root.parent).free>=plan['resources']['minimum_free_disk_gib']*2**30
    root.mkdir(exist_ok=True);lock=(root/'stage.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (root/'receipt.json').exists(),'Completed full fitting cannot restart'
    state=dict(plan_sha256=sha(path),source_contract=source['fit_contract'],schema=SCHEMA)
    marker=root/'stage_plan.json'
    if marker.exists():assert json.loads(marker.read_text())==state
    else:marker.write_text(json.dumps(state,indent=2)+'\n')
    folder=root/'cohorts';folder.mkdir(exist_ok=True);manifest=[];counts=Counter();total=0
    for number,(cohort,rows,entries) in enumerate(cohorts(source,plan),1):
        stem=cohort['cohort_id'];part=folder/(stem+'.jsonl.gz');cp=folder/(stem+'.receipt.json');started=time.perf_counter()
        if cp.exists():
            saved=json.loads(cp.read_text());assert saved['stage']==state and saved['cohort_id']==stem
            assert saved['candidate_sha256']==sha(part)
            with gzip.open(part,'rt') as f:
                exported=(json.loads(line) for line in f)
                local=Counter();seen=0
                for expected,_,_,_,_,_ in cases(source,plan,cohort,entries):
                    value=next(exported)
                    assert all(value[k]==v for k,v in expected.items())
                    local[value['disposition']]+=1;seen+=1
                assert next(exported,None) is None
            assert saved['candidate_rows']==seen and saved['status_counts']==dict(local)
        else:
            operators={}
            for _,_,_,qualified in entries:
                for (mode,tree),(audit,_,_) in qualified.items():
                    if mode not in operators:operators[mode]=operators_for(source,rows,audit)
            temporary=part.with_suffix('.gz.partial');local=Counter();seen=0
            with gzip.open(temporary,'wt') as output:
                for expected,matrix,y,audit,original,cert in cases(source,plan,cohort,entries):
                    value=candidate(source,plan,rows,expected,matrix,y,operators[expected['loading_mode']],audit,original,cert)
                    output.write(json.dumps(value,sort_keys=True,allow_nan=False)+'\n');local[value['disposition']]+=1;seen+=1
            with temporary.open('rb') as f:os.fsync(f.fileno())
            os.replace(temporary,part)
            saved=dict(stage=state,cohort_id=stem,candidate_rows=seen,status_counts=dict(local),candidate_sha256=sha(part))
            temporary_receipt=cp.with_suffix('.json.partial');temporary_receipt.write_text(json.dumps(saved,indent=2)+'\n')
            with temporary_receipt.open('rb') as f:os.fsync(f.fileno())
            os.replace(temporary_receipt,cp)
        counts.update(local);total+=seen
        manifest.append(dict(cohort_id=stem,path=str(part.relative_to(root)),sha256=sha(part),
            receipt_path=str(cp.relative_to(root)),receipt_sha256=sha(cp),candidate_rows=seen))
        print('full_shared_entity_fit_cohort',number,'/',len(source['cohorts']),'records',len(rows),
            'candidates',seen,'seconds',time.perf_counter()-started,'checkpoint',cp.exists(),flush=True)
        if stop_after_cohorts==number:raise InterruptedError('Software fitting checkpoint contract')
    (root/'cohort_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    database=root/'candidate_index.sqlite'
    if database.exists():database.unlink()  # This index is rebuilt from immutable cohort receipts.
    db=sqlite3.connect(database);db.execute('PRAGMA cache_size=-65536')
    db.execute('CREATE TABLE candidates(id TEXT PRIMARY KEY,fit TEXT,mode TEXT,tree TEXT,method TEXT,audit TEXT,status TEXT,UNIQUE(fit,mode,tree,method))')
    for part in manifest:
        with gzip.open(root/part['path'],'rt') as f:
            values=[json.loads(line) for line in f]
        db.executemany('INSERT INTO candidates VALUES (?,?,?,?,?,?,?)',[(v['candidate_id'],v['fit_input_id'],v['loading_mode'],v['tree'],v['method'],v['covariance_audit_id'],v['disposition']) for v in values])
    db.commit();links=0;settings=Counter();fields=SETTING_FIELDS+NEW_LINK_EXTRA+FIT_LINK_EXTRA
    with gzip.open(source['qualification_root']/'setting_audit_links.tsv.gz','rt') as f,gzip.open(root/'setting_fit_links.tsv.gz','wt') as g:
        reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==SETTING_FIELDS+NEW_LINK_EXTRA
        writer=csv.DictWriter(g,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in reader:
            for method in plan['methods']:
                cid=candidate_id(source['fit_contract'],row['fit_input_id'],row['loading_mode'],row['tree'],method)
                found=db.execute('SELECT audit,status FROM candidates WHERE id=?',(cid,)).fetchone()
                assert found is not None and found[0]==row['audit_id']
                writer.writerow({**row,'likelihood_method':method,'candidate_id':cid,'candidate_disposition':found[1]})
                links+=1;settings[found[1]]+=1
    db.close();database.unlink()
    q=source['qualification_completion'];unique_inputs=q['unique_designs']*2
    assert total==unique_inputs*len(plan['methods'])*len(plan['loading_modes'])*len(plan['trees'])
    assert links==q['setting_audit_links']*len(plan['methods'])
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=q['model_setting_rows'],unique_cohorts=len(manifest),
        unique_designs=q['unique_designs'],unique_fit_inputs=unique_inputs,candidate_rows=total,setting_fit_links=links,
        candidate_status_counts=dict(counts),setting_status_counts=dict(settings),methods=plan['methods'],trees=plan['trees'],loading_modes=plan['loading_modes'])
    verify(bindings)
    names=['stage_plan.json','cohort_manifest.json','setting_fit_links.tsv.gz']+[v[k] for v in manifest for k in ['path','receipt_path']]
    receipt=dict(status=PRODUCER_STATUS,plan_sha256=sha(path),
        source_contract=source['fit_contract'],**summary,artifacts={name:sha(root/name) for name in names},
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with (root/'receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
