#!/usr/bin/env python3
"""Reconstruct every case and selection membership with an independent SQL census."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile
from full_matching_case_sources import load, MASKS, BIT_FIELDS, CASE_FIELDS, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


class BitUnion:
    def __init__(self): self.value=0
    def step(self,value): self.value |= value
    def finalize(self): return self.value


def digest(tag,a,b):
    return hashlib.sha256(json.dumps([tag,a,b],separators=(',',':')).encode()).hexdigest()


def run(path, output):
    path=Path(path);plan=json.loads(path.read_text());source,bindings=load(plan,path);source_bindings=dict(bindings)
    root=Path(plan['output']);rp=root/'receipt.json';r=json.loads(rp.read_text())
    assert r['status']=='complete_full_matching_case_index_pending_independent_readback' and r['plan_sha256']==sha(path) and r['scientific_eligibility'] is False
    assert r['source_hashes']==source_bindings
    unmatched_path=source['root']/'target_policy_coverage_status.tsv.gz'
    assert r['original_unmatched_status_path']==str(unmatched_path) and r['original_unmatched_status_sha256']==sha(unmatched_path)
    bind(bindings,rp)
    assert set(r['artifacts'])=={'stage_plan.json','selection_case_links.tsv.gz','case_index.tsv.gz'}
    for n,h in r['artifacts'].items():bind(bindings,root/n,h)
    verify(bindings)
    assert json.loads((root/'stage_plan.json').read_text())==dict(plan_sha256=sha(path),schema='fixed-matched-logical-case-v1',recovery='Full deterministic replay; no partial census is accepted.')
    policies={p:i for i,p in enumerate(plan['policies'])};scenarios={s['scenario_id']:i for i,s in enumerate(source['scenarios'])}
    with tempfile.TemporaryDirectory(prefix='independent-full-case-index-',dir=root) as directory:
        db=sqlite3.connect(Path(directory)/'source.sqlite');db.create_aggregate('bit_union',1,BitUnion)
        db.execute('PRAGMA temp_store=FILE');db.execute('PRAGMA cache_size=-262144')
        db.execute('CREATE TABLE policies(tid TEXT,pidx INTEGER,guide TEXT,matched INTEGER,full INTEGER,plddt70 INTEGER,both INTEGER,PRIMARY KEY(tid,pidx))')
        db.execute('CREATE TABLE chosen(tid TEXT,bid TEXT,pidx INTEGER,sidx INTEGER,endpoint INTEGER,guide TEXT,physical TEXT,tf INTEGER,bf INTEGER,jf INTEGER,tp INTEGER,bp INTEGER,jp INTEGER,tb INTEGER,bb INTEGER,jb INTEGER,PRIMARY KEY(tid,pidx,sidx))')
        statuses=unmatched=expected_selections=0;membership={};targetflags={}
        with gzip.open(unmatched_path,'rt') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                tid=row['target_id'];pidx=policies[row['policy']];node=source['nodes']['target'][tid]
                assert row['guide']==node['guide'] and row['family']==node['family'] and row['focal_taxon']==node['taxon_id'] and row['gene_node']==node['gene_node'] and row['target_pair_key']==node['pair_key']
                have=row['matched_scenarios'].split(',') if row['matched_scenarios'] else []
                missing=row['unmatched_scenarios'].split(',') if row['unmatched_scenarios'] else []
                assert len(have)==len(set(have)) and len(missing)==len(set(missing)) and set(have).isdisjoint(missing) and set(have+missing)==set(scenarios)
                flag=0
                for sid in have:flag |= 2**scenarios[sid]
                key=tid,pidx;assert key not in membership;membership[key]=flag
                bits=tuple(int(row[f'target_{m}_pass_bits']) for m in MASKS)
                assert all(0<=v<64 for v in bits) and bits[2]==bits[0]&bits[1]
                if node['same_model']:assert bits==(0,0,0)
                assert tid not in targetflags or targetflags[tid]==bits;targetflags[tid]=bits
                db.execute('INSERT INTO policies VALUES(?,?,?,?,?,?,?)',(tid,pidx,node['guide'],flag,*bits));statuses+=1;unmatched+=len(missing);expected_selections+=len(have)
        assert len(membership)==plan['expected']['target_policy_records']==len(source['nodes']['target'])*len(policies)
        db.commit();selected=retained=0;pending=[]
        with gzip.open(source['root']/'selected_pair_coverage.tsv.gz','rt') as a,gzip.open(root/'selection_case_links.tsv.gz','rt') as b:
            original=csv.DictReader(a,delimiter='\t');exported=csv.DictReader(b,delimiter='\t')
            assert exported.fieldnames==['source_row_ordinal','case_id','physical_case_id']+original.fieldnames
            for raw,actual in itertools.zip_longest(original,exported):
                assert raw is not None and actual is not None
                tid,bid=raw['target_id'],raw['background_id'];t=source['nodes']['target'][tid];bg=source['nodes']['background'][bid]
                pidx,sidx=policies[raw['policy']],scenarios[raw['scenario_id']]
                assert membership[tid,pidx] & 2**sidx and raw['endpoint_order'] in ['0','1']
                assert t['guide']==bg['guide']==raw['guide'] and raw['family']==t['family'] and raw['focal_taxon']==t['taxon_id'] and raw['gene_node']==t['gene_node']
                assert raw['target_pair_key']==t['pair_key'] and raw['background_pair_key']==bg['pair_key']
                assert raw['target_sequence_distance']==str(t['sequence_distance']) and raw['background_sequence_distance']==str(bg['sequence_distance'])
                for n in [t,bg]:
                    ends=sorted([(n['model_id_'+s],n['version_'+s]) for s in ['a','b']])
                    assert hashlib.sha256(json.dumps(ends,separators=(',',':')).encode()).hexdigest()==n['pair_key'] and int(ends[0]==ends[1])==n['same_model']
                numbers=[int(raw[f]) for f in BIT_FIELDS]
                for i,m in enumerate(MASKS):
                    tf,bf,jf=numbers[3*i:3*i+3];assert tf==targetflags[tid][i] and 0<=bf<64 and jf==tf&bf
                    retained+=jf.bit_count()
                for i in range(3):assert numbers[6+i]==numbers[i]&numbers[3+i]
                if bg['same_model']:assert all(numbers[i]==0 for i in [1,4,7])
                selected+=1;cid=digest('fixed-matched-logical-case-v1',tid,bid);pid=digest('fixed-matched-physical-case-v1',t['pair_key'],bg['pair_key'])
                assert actual==dict(source_row_ordinal=str(selected),case_id=cid,physical_case_id=pid,**raw)
                pending.append((tid,bid,pidx,sidx,int(raw['endpoint_order']),t['guide'],pid,*numbers))
                if len(pending)==10000:db.executemany('INSERT INTO chosen VALUES('+','.join(['?']*16)+')',pending);pending=[]
                if selected%250000==0:print('independent_full_matching_case_selections',selected,flush=True)
        if pending:db.executemany('INSERT INTO chosen VALUES('+','.join(['?']*16)+')',pending)
        db.commit();assert selected==expected_selections==plan['expected']['selected_records'] and unmatched==plan['expected']['unmatched_decisions']
        # Known unique policy/scenario keys, legal membership and exact total imply
        # every declared selected scenario is present; no vanished choice can hide.
        db.create_function('bit_count',1,int.bit_count)
        assert not db.execute('SELECT p.tid,p.pidx FROM policies p LEFT JOIN chosen c ON p.tid=c.tid AND p.pidx=c.pidx GROUP BY p.tid,p.pidx HAVING COUNT(c.sidx)<>bit_count(p.matched)').fetchall()
        columns=['tf','bf','jf','tp','bp','jp','tb','bb','jb']
        query='SELECT tid,bid,COUNT(*),bit_union(1<<endpoint),bit_union(1<<pidx),bit_union(1<<sidx),'+','.join('MIN('+c+'),MAX('+c+')' for c in columns)+' FROM chosen GROUP BY tid,bid ORDER BY tid,bid'
        disposition=Counter();physical=set();ts=set();bs=set();logical=0;maximum=0
        guide_sets={g:{k:set() for k in ['physical_cases','target_nodes','background_nodes','target_families','background_families','focal_taxa']} for g in plan['guides']}
        guide_counts={g:dict(logical_cases=0,selected_records=0) for g in plan['guides']}
        with gzip.open(root/'case_index.tsv.gz','rt') as handle:
            table=csv.DictReader(handle,delimiter='\t');assert table.fieldnames==CASE_FIELDS
            for sql,actual in itertools.zip_longest(db.execute(query),table):
                assert sql is not None and actual is not None
                tid,bid,n,eb,pb,sb=sql[:6];t,bg=source['nodes']['target'][tid],source['nodes']['background'][bid]
                values=sql[6:];assert all(values[2*i]==values[2*i+1] for i in range(9));bits=values[::2]
                pid=digest('fixed-matched-physical-case-v1',t['pair_key'],bg['pair_key'])
                expected=dict(case_id=digest('fixed-matched-logical-case-v1',tid,bid),physical_case_id=pid,target_id=tid,background_id=bid,guide=t['guide'],target_family=t['family'],background_family=bg['family'],focal_taxon=t['taxon_id'],gene_node=t['gene_node'],target_gene_a=t['gene_a'],target_gene_b=t['gene_b'],background_gene_a=bg['gene_a'],background_gene_b=bg['gene_b'],background_taxon_a=bg['taxon_a'],background_taxon_b=bg['taxon_b'],target_pair_key=t['pair_key'],background_pair_key=bg['pair_key'],target_same_model=t['same_model'],background_same_model=bg['same_model'],target_sequence_distance=t['sequence_distance'],background_sequence_distance=bg['sequence_distance'],selection_records=n,endpoint_order_bits=eb,policy_bits=pb,scenario_bits=sb,**dict(zip(BIT_FIELDS,bits)))
                assert actual=={k:str(v) for k,v in expected.items()}
                logical+=1;maximum=max(maximum,n);physical.add(pid);ts.add(tid);bs.add(bid)
                disposition[f"target_same_model={t['same_model']},background_same_model={bg['same_model']}"]+=1
                g=t['guide'];guide_counts[g]['logical_cases']+=1;guide_counts[g]['selected_records']+=n
                for k,v in [('physical_cases',pid),('target_nodes',tid),('background_nodes',bid),('target_families',t['family']),('background_families',bg['family']),('focal_taxa',t['taxon_id'])]:guide_sets[g][k].add(v)
        guides={g:dict(**guide_counts[g],**{k:len(v) for k,v in guide_sets[g].items()}) for g in plan['guides']}
        assert retained==sum(int(row['joint_pass_matched_records']) for row in source['attrition'].values())
        summary=dict(target_nodes=len(source['nodes']['target']),background_nodes=len(source['nodes']['background']),target_policy_records=statuses,scenarios=len(scenarios),scenario_decisions=statuses*len(scenarios),selected_records=selected,unmatched_decisions=unmatched,logical_cases=logical,physical_cases=len(physical),unique_selected_targets=len(ts),unique_selected_backgrounds=len(bs),maximum_selection_reuse=maximum,case_dispositions=dict(disposition),guide_census=guides,future_case_mask_rows=2*logical,future_order_pair_cells=8*logical,retained_selection_screen_cells=retained,screens=plan['screens'],masks=MASKS,guides=plan['guides'],policies=plan['policies'])
        assert all(summary[k]==r[k] for k in SUMMARY_FIELDS)
        db.close()
    verify(bindings)
    result=dict(status='passed_full_matching_case_index_sql_readback',plan_sha256=sha(path),producer_receipt_sha256=sha(rp),**summary,source_hashes=bindings,scientific_eligibility=False,scope=plan['scope'])
    with Path(output).open('x') as handle:handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.plan,a.output)
