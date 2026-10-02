#!/usr/bin/env python3
"""Run both sequence aligners and every input order on all closed physical sets."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import json
import fcntl
from pathlib import Path
import shutil
import sqlite3
import zlib

from full_triad_sequence_sources import load_catalog
from reference_measurement_union_sources import verify
from run_ortholog_pair_guide_comparison import sha
from triad_sequence_native import METHODS, ORDERS, digest, execute, validate_checkpoint


def run(plan_path):
    plan=json.loads(Path(plan_path).read_text());ph=sha(plan_path)
    assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    sets,source,bindings=load_catalog(plan,plan_path);index={r['sequence_set_id']:r for r in sets}
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    configuration=out/'configuration.json'
    if configuration.exists():assert json.loads(configuration.read_text())['plan_sha256']==ph
    else:
        assert set(out.iterdir())=={out/'run.lock'}
        configuration.write_text(json.dumps(dict(plan_sha256=ph,methods=METHODS,orders=ORDERS),indent=2)+'\n')
    assert not (out/'receipt.json').exists(),'Completed run cannot be restarted'
    db=sqlite3.connect(out/'native_alignments.sqlite')
    db.execute('CREATE TABLE IF NOT EXISTS alignments (sequence_set_id TEXT, method TEXT, permutation TEXT, payload BLOB NOT NULL, payload_sha256 TEXT NOT NULL, PRIMARY KEY(sequence_set_id,method,permutation))')
    existing=set();counts=Counter();rows=0;common=0
    for sid,method,permutation,payload,h in db.execute('SELECT * FROM alignments'):
        raw=zlib.decompress(payload);assert digest(raw)==h;result=json.loads(raw)
        assert (sid,method,permutation)==tuple(result[k] for k in ['sequence_set_id','method','permutation'])
        validate_checkpoint(plan,index[sid],result,out)
        existing.add((sid,method,permutation));counts[method+':'+result['status']]+=1;rows+=1
        common+=result['common_residues'] or 0
    def jobs():
        for record in sets:
            for method in METHODS:
                for order in ORDERS:
                    key=record['sequence_set_id'],method,''.join(map(str,order))
                    if key not in existing:yield record,method,order
    iterator=iter(jobs());pending={}
    def refill(pool):
        while len(pending)<2*plan['workers']:
            job=next(iterator,None)
            if job is None:break
            pending[pool.submit(execute,plan,*job,out)]=job
    def checkpoint():
        db.commit()
        state=dict(status='running_full_triad_sequence_alignment_dispositions',completed=rows,
            total=12*len(sets),counts=dict(counts),plan_sha256=ph)
        temp=out/'state.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(out/'state.json')
        assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
        assert (out/'native_alignments.sqlite').stat().st_size<=plan['resources']['output_allowance_gib']*2**30
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        refill(pool)
        while pending:
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                record,method,order=pending.pop(future);result=future.result()
                validate_checkpoint(plan,record,result,out)
                raw=json.dumps(result,separators=(',',':')).encode()
                db.execute('INSERT INTO alignments VALUES (?,?,?,?,?)',
                    (record['sequence_set_id'],method,result['permutation'],zlib.compress(raw,1),digest(raw)))
                rows+=1;counts[method+':'+result['status']]+=1;common+=result['common_residues'] or 0
                if rows%128==0:
                    checkpoint();print('Full sequence-control native states',rows,'/',12*len(sets),flush=True)
            refill(pool)
    checkpoint();assert rows==source['native_alignment_states']==12*len(sets)
    db.close();verify(bindings)
    result=dict(status='complete_full_triad_sequence_alignment_dispositions_pending_independent_readback',
        plan_sha256=ph,ordered_model_triads=source['ordered_model_triads'],source_ready_triads=source['source_ready_triads'],
        unique_sequence_model_sets=len(sets),native_alignment_states=rows,counts=dict(counts),
        common_residue_occurrences=common,source_hashes=bindings,
        artifacts={name:sha(out/name) for name in ['configuration.json','native_alignments.sqlite','state.json']},
        scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
