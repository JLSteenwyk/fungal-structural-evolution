#!/usr/bin/env python3
"""Distinguish aligned-core, complete-interval and full-chain sequence identity."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import Counter


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def textsha(s):return hashlib.sha256(s.encode()).hexdigest()
def rows(p):
    with Path(p).open() as f:yield from csv.DictReader(f,delimiter='\t')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
    for f,h in p['pins'].items():assert sha(f)==h
    out=Path(p['output']);assert not out.exists()
    candidates=list(rows(p['candidates']));needed={t for x in candidates for t in json.loads(x['triad_keys_json'])};triads={}
    for line in Path(p['triads']).open():
        r=json.loads(line)
        if r['triad_key'] in needed:triads[r['triad_key']]=r
    assert set(triads)==needed
    ids={x[k] for x in triads.values() for k in ['interval_a','interval_b']};intervals={};inputs={}
    for line in Path(p['intervals']).open():
        x=json.loads(line)
        if x['interval_id'] in ids:intervals[x['interval_id']]=x
    for line in Path(p['inputs']).open():
        x=json.loads(line)
        if x['interval_id'] in ids and x['mask']=='full':inputs[x['interval_id']]=x
    assert set(intervals)==set(inputs)==ids
    for iid in ids:
        i=intervals[iid];x=inputs[iid];assert x['status']=='ready'
        assert all(i[k]==x[k] for k in ['model_id','version','start','end','original_length','source_sha256'])
        assert len(x['sequence'])==i['length']==i['end']-i['start']+1
    details=[];bytriad={}
    for tk,t in sorted(triads.items()):
        r=dict(triad_key=tk)
        for role in ['a','b']:
            iid=t['interval_'+role];i=intervals[iid];x=inputs[iid]
            r['interval_'+role]=iid
            for k in ['model_id','version','start','end','length','original_length','sequence_sha256','source_sha256']:r[role+'_'+k]=i[k]
            r[role+'_interval_sequence_sha256']=textsha(x['sequence'])
        r['complete_intervals_identical']=int(inputs[t['interval_a']]['sequence']==inputs[t['interval_b']]['sequence'])
        r['full_chain_sequence_hashes_identical']=int(r['a_sequence_sha256']==r['b_sequence_sha256'])
        r['same_model_version']=int((r['a_model_id'],r['a_version'])==(r['b_model_id'],r['b_version']))
        details.append(r);bytriad[tk]=r
    output=[]
    for x in candidates:
        pool=[bytriad[t] for t in json.loads(x['triad_keys_json'])];n=len(pool)
        models={(r['a_model_id'],r['a_version'],r['b_model_id'],r['b_version']) for r in pool};assert len(models)==1
        identical_core=float(x['all_alternatives_min_sequence_identity_ab'])==1
        domain=sum(r['complete_intervals_identical'] for r in pool);chain={r['full_chain_sequence_hashes_identical'] for r in pool};same={r['same_model_version'] for r in pool};assert len(chain)==len(same)==1
        r=dict(x,all_aligned_cores_identical=int(identical_core),complete_interval_pairs=n,identical_complete_interval_pairs=domain,all_complete_intervals_identical=int(domain==n),full_chain_sequence_hashes_identical=next(iter(chain)),same_model_version=next(iter(same)))
        for role in ['a','b']:
            for field in ['model_id','version','original_length','sequence_sha256']:r[role+'_'+field]=pool[0][role+'_'+field]
        output.append(r)
    out.mkdir(parents=True)
    for name,data in [('candidate_sequence_context.tsv',output),('interval_pair_context.tsv',details)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    counts=Counter((x['candidate_class'],x['study_role'],x['all_aligned_cores_identical'],x['all_complete_intervals_identical'],x['full_chain_sequence_hashes_identical']) for x in output)
    receipt=dict(status='complete_candidate_sequence_context_pending_readback',plan_sha256=sha(a.plan),candidates=len(output),triads=len(details),intervals=len(ids),strata=[dict(candidate_class=k[0],study_role=k[1],all_aligned_cores_identical=k[2],all_complete_intervals_identical=k[3],full_chain_sequence_hashes_identical=k[4],count=v) for k,v in sorted(counts.items())],artifacts={f.name:sha(f) for f in out.iterdir()},scope='Aligned-core identity from audited common-core fits; complete-interval exact string equality from full-mask inputs; full-chain identity from source sequence hashes. Different intervals can reflect boundary/context differences. These are not causal attribution, accuracy or sequence-independent evolutionary evidence.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
