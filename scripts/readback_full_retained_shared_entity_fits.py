#!/usr/bin/env python3
"""Full original-fit replay with component-local independent spectral audits."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import json
from pathlib import Path
import sqlite3
from full_retained_shared_entity_fit_sources import load,cohorts,cases,candidate_id,FIT_LINK_EXTRA,SUMMARY_FIELDS,NEW_LINK_EXTRA,SCHEMA,PRODUCER_STATUS,READER_STATUS,operators_for
from readback_retained_shared_entity_candidate import numeric
from full_expanded_model_design_sources import SETTING_FIELDS,digest
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def run(path,output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);original=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';receipt=json.loads(rp.read_text())
    lock=(root/'readback.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    completed_marker=root/'independent_readback_completed.json'
    if completed_marker.exists():
        previous=json.loads(completed_marker.read_text());assert sha(previous['output'])==previous['output_sha256']
        raise AssertionError('Completed independent fitting readback cannot rerun under another output name')
    assert not Path(output).exists(),'Completed independent fitting readback cannot restart'
    assert receipt['status']==PRODUCER_STATUS
    assert receipt['plan_sha256']==sha(path) and receipt['source_contract']==source['fit_contract']
    assert receipt['source_hashes']==original and receipt['scientific_eligibility'] is False
    bind(bindings,rp)
    for name,d in receipt['artifacts'].items():bind(bindings,root/name,d)
    verify(bindings)
    stage=dict(plan_sha256=sha(path),source_contract=source['fit_contract'],schema=SCHEMA)
    assert json.loads((root/'stage_plan.json').read_text())==stage
    readback_state=dict(schema='full-retained-shared-entity-component-spectral-readback-v1',plan_sha256=sha(path),producer_receipt_sha256=sha(rp))
    marker=root/'readback_stage.json'
    if marker.exists():assert json.loads(marker.read_text())==readback_state
    else:marker.write_text(json.dumps(readback_state,indent=2)+'\n')
    manifest=json.loads((root/'cohort_manifest.json').read_text())
    assert [r['cohort_id'] for r in manifest]==[r['cohort_id'] for r in source['cohorts']]
    assert len({r['cohort_id'] for r in manifest})==len(manifest)
    expected_artifacts={'stage_plan.json','cohort_manifest.json','setting_fit_links.tsv.gz'}|{r[k] for r in manifest for k in ['path','receipt_path']}
    assert set(receipt['artifacts'])==expected_artifacts
    database=root/'fit_readback.sqlite'
    if database.exists():database.unlink()  # Derived scratch; lock/state guard a full independent replay.
    db=sqlite3.connect(database);db.execute('PRAGMA cache_size=-65536')
    db.execute('CREATE TABLE candidates(id TEXT PRIMARY KEY,fit TEXT,mode TEXT,tree TEXT,method TEXT,audit TEXT,status TEXT,independent TEXT,UNIQUE(fit,mode,tree,method))')
    counts=Counter();audited=Counter();total=0;design_count=0
    audits_path=root/'independent_candidate_audits.jsonl.gz'
    with gzip.open(audits_path,'wt') as audit_output:
        for number,((cohort,rows,entries),part) in enumerate(zip(cohorts(source,plan),manifest),1):
            cid=cohort['cohort_id'];assert part['path']=='cohorts/'+cid+'.jsonl.gz' and part['receipt_path']=='cohorts/'+cid+'.receipt.json'
            saved=json.loads((root/part['receipt_path']).read_text());assert saved['stage']==stage and saved['cohort_id']==cid
            assert saved['candidate_sha256']==part['sha256']==receipt['artifacts'][part['path']]
            assert part['receipt_sha256']==receipt['artifacts'][part['receipt_path']]
            operators={}
            for _,_,_,qualified in entries:
                for (mode,tree),(audit,_,_) in qualified.items():
                    if mode not in operators:operators[mode]=operators_for(source,rows,audit)
            local=Counter();seen=0
            with gzip.open(root/part['path'],'rt') as f:
                exported=(json.loads(line) for line in f)
                for expected,matrix,response,source_audit,original,cert in cases(source,plan,cohort,entries):
                    value=next(exported);assert value['scientific_eligibility'] is False
                    audit=numeric(source,plan,rows,expected,matrix,response,value,operators.get(expected['loading_mode']),
                        source_audit,original,cert)
                    record=dict(candidate_id=expected['candidate_id'],source_candidate_sha256=digest(value),**audit)
                    audit_output.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
                    db.execute('INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?)',(value['candidate_id'],value['fit_input_id'],value['loading_mode'],value['tree'],value['method'],value['covariance_audit_id'],value['disposition'],audit['disposition']))
                    local[value['disposition']]+=1;audited[audit['disposition']]+=1;total+=1;seen+=1
                assert next(exported,None) is None
            assert part['candidate_rows']==saved['candidate_rows']==seen and saved['status_counts']==dict(local)
            counts.update(local);design_count+=len(entries);db.commit()
            print('independent_full_shared_entity_cohort',number,'/',len(manifest),'candidates',seen,flush=True)
    settings=Counter();links=0;link_audits=Counter()
    independent_links=root/'independent_setting_fit_links.tsv.gz'
    with gzip.open(source['qualification_root']/'setting_audit_links.tsv.gz','rt') as f,gzip.open(root/'setting_fit_links.tsv.gz','rt') as g,gzip.open(independent_links,'wt') as h:
        original_rows=csv.DictReader(f,delimiter='\t');reader=csv.DictReader(g,delimiter='\t')
        assert original_rows.fieldnames==SETTING_FIELDS+NEW_LINK_EXTRA and reader.fieldnames==SETTING_FIELDS+NEW_LINK_EXTRA+FIT_LINK_EXTRA
        writer=csv.DictWriter(h,reader.fieldnames+['independent_disposition'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in original_rows:
            for method in plan['methods']:
                exported=next(reader);cid=candidate_id(source['fit_contract'],row['fit_input_id'],row['loading_mode'],row['tree'],method)
                saved=db.execute('SELECT audit,status,independent FROM candidates WHERE id=?',(cid,)).fetchone()
                assert saved is not None and saved[0]==row['audit_id']
                assert exported=={**row,'likelihood_method':method,'candidate_id':cid,'candidate_disposition':saved[1]}
                writer.writerow({**exported,'independent_disposition':saved[2]})
                settings[saved[1]]+=1;link_audits[saved[2]]+=1;links+=1
        assert next(reader,None) is None
    db.close();database.unlink()
    q=source['qualification_completion']
    assert total==design_count*2*len(plan['methods'])*len(plan['loading_modes'])*len(plan['trees'])
    assert design_count==q['unique_designs'] and links==q['setting_audit_links']*len(plan['methods'])
    summary=dict(logical_cases=len(source['ids']),model_setting_rows=q['model_setting_rows'],unique_cohorts=len(manifest),
        unique_designs=design_count,unique_fit_inputs=design_count*2,candidate_rows=total,setting_fit_links=links,
        candidate_status_counts=dict(counts),setting_status_counts=dict(settings),methods=plan['methods'],trees=plan['trees'],loading_modes=plan['loading_modes'])
    assert all(receipt[k]==summary[k] for k in SUMMARY_FIELDS)
    for p in [audits_path,independent_links,marker]:bind(bindings,p)
    verify(bindings)
    result=dict(status=READER_STATUS,plan_sha256=sha(path),producer_receipt_sha256=sha(rp),
        source_contract=source['fit_contract'],independent_backend=plan['independent_backend'],**summary,independent_candidate_status_counts=dict(audited),
        independent_setting_status_counts=dict(link_audits),source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    with completed_marker.open('x') as f:json.dump(dict(output=str(output),output_sha256=sha(output),
        plan_sha256=sha(path),producer_receipt_sha256=sha(rp)),f,indent=2);f.write('\n')
    print(json.dumps(summary),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
