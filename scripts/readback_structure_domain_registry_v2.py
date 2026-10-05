#!/usr/bin/env python3
"""Independently reconstruct all registry rows from frozen source records."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import json
from pathlib import Path
import sqlite3
import time
import psutil
from catalog_whole_proteome_structures import sha


def expected_rows(payload, agree):
    segments={};members=set()
    for policy,entry in payload['policies'].items():
        alternative=payload['alternatives'][entry['alternative_index']]
        annotations=alternative['annotations']
        spans=[(h['alignment_start'],h['alignment_end']) for h in annotations]
        overlap=any(max(a,c)<=min(b,d) for i,(a,b) in enumerate(spans) for c,d in spans[i+1:])
        whole_ok=bool(agree) and not overlap and all(Decimal(h['hmm_coverage'])>=Decimal('0.70') for h in annotations) and all(entry[k]==0 for k in ['primary_rank_tie_pairs','unresolved_overlap_pairs','candidate_nested_pairs'])
        for h in annotations:
            hit=h['hit_id'];fields=(hit,h['pfam_accession'],h['pfam_type'],h['pfam_clan'],h['alignment_start'],h['alignment_end'],h['envelope_start'],h['envelope_end'],str(Decimal(h['hmm_coverage'])))
            if hit in segments and segments[hit]!=fields:raise ValueError('Inconsistent source hit')
            segments[hit]=fields
            candidate=whole_ok and h['pfam_type']=='Domain' and Decimal(h['hmm_coverage'])>=Decimal('0.70') and h['alignment_end']-h['alignment_start']+1>=30
            record=(policy,hit,int(whole_ok),int(candidate))
            if record in members:raise ValueError('Duplicate source policy hit')
            members.add(record)
    return set(segments.values()),members


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--receipt',type=Path,required=True)
    args=ap.parse_args();assert not args.receipt.exists();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Readback plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Readback dependency changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    def state(stage,**kw):(out/'state.json').write_text(json.dumps(dict(stage=stage,**kw))+'\n')
    state('verifying_completed_original_registry_transport')
    adapter_path=Path(plan['producer_receipt'])
    transport_path=Path(plan['producer_transport'])
    adapter,transport=[json.loads(path.read_text()) for path in [adapter_path,transport_path]]
    assert adapter['status']=='completed_full_refreshed_afdb_domain_registry_pending_independent_readback'
    assert transport['validation_sha256']==sha(adapter_path)
    assert transport['original_tool_terminal_exit_code']==0
    assert transport['whole_wrapper_initial_and_terminal_payloads_matched']
    assert transport['manager_start_records']==transport['manager_completion_records']==1
    for mapping in [adapter['source_hashes'],transport['source_hashes']]:
        for path,digest in mapping.items():
            assert sha(path)==digest,path
    pins={str(args.plan):ph,**plan['pins'],**adapter['source_hashes']}
    for path in [adapter_path,transport_path,Path(__file__)]:pins[str(path)]=sha(path)
    verify();source=json.loads(Path(plan['producer_plan']).read_text());root=Path(source['output']);receipt=json.loads((root/'receipt.json').read_text())
    if receipt['status']!='complete_structure_domain_interval_registry_pending_independent_readback' or receipt['plan_sha256']!=sha(plan['producer_plan']):raise ValueError('Registry not complete or wrong plan')
    database=root/'structure_domains.sqlite'
    if sha(database)!=receipt['database_sha256']:raise ValueError('Registry hash differs')
    db=sqlite3.connect('file:'+str(database.resolve())+'?mode=ro',uri=True)
    ann=sqlite3.connect('file:'+str(Path(source['architectures']).resolve())+'?mode=ro',uri=True)
    counts=Counter();candidates=Counter();seen=set();started=time.monotonic()
    with Path(source['models']).open() as f:
        for line in f:
            model=json.loads(line);key=Path(model['path']).stem
            if key in seen:raise ValueError('Duplicate source model')
            seen.add(key)
            raw,agree,text=ann.execute('SELECT raw_hits,retained_sets_agree,candidate_architectures_json FROM queries WHERE sequence_id=?',('S'+model['sequence_sha256'],)).fetchone()
            actual=db.execute('SELECT * FROM models WHERE model_key=?',(key,)).fetchone()
            expected=(key,model['model_id'],model['version'],model['sequence_sha256'],model['length'],model['path'],raw,agree)
            if actual!=expected:raise ValueError('Model metadata differs: '+key)
            segments,members=expected_rows(json.loads(text),agree)
            actual_segments=set(db.execute('SELECT hit_id,pfam_accession,pfam_type,pfam_clan,alignment_start,alignment_end,envelope_start,envelope_end,hmm_coverage FROM segments WHERE model_key=?',(key,)))
            actual_members=set(db.execute('SELECT policy,hit_id,conservative_architecture,domain_interval_candidate FROM policy_hits WHERE model_key=?',(key,)))
            if segments!=actual_segments or members!=actual_members:raise ValueError('Interval or policy rows differ: '+key)
            counts['models']+=1;counts['models_without_hits']+=raw==0;counts['segments']+=len(segments);counts['policy_hits']+=len(members)
            for policy,_,_,eligible in members:candidates[policy]+=eligible
            if counts['models']%10000==0:state('checking_all_models',counts=dict(counts),elapsed_seconds=time.monotonic()-started)
    seen_proteins=set()
    with Path(source['protein_links']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            identity=(row['taxon_id'],row['protein_id']);key=Path(row['model_path']).stem
            if identity in seen_proteins:raise ValueError('Duplicate source protein')
            seen_proteins.add(identity)
            actual=db.execute('SELECT p.model_key,m.sequence_sha256,m.model_id,m.version FROM protein_links p JOIN models m USING(model_key) WHERE p.taxon_id=? AND p.protein_id=?',identity).fetchone()
            if actual!=(key,row['sequence_sha256'],row['model_id'],int(row['version'])):raise ValueError('Protein/model link differs')
            if ann.execute('SELECT sequence_id FROM proteins WHERE taxon_id=? AND protein_id=?',identity).fetchone()!=('S'+row['sequence_sha256'],):raise ValueError('Annotation protein identity differs')
            counts['protein_links']+=1
    for table in ['models','segments','policy_hits','protein_links']:
        if db.execute('SELECT count(*) FROM '+table).fetchone()[0]!=counts[table]:raise ValueError('Extra or missing database rows')
    if dict(counts)!=receipt['counts'] or dict(candidates)!=receipt['candidate_domain_intervals_by_policy']:raise ValueError('Receipt totals differ')
    if db.execute('PRAGMA integrity_check').fetchone()!=('ok',) or db.execute('PRAGMA foreign_key_check').fetchall():raise ValueError('Registry integrity failure')
    db.close();ann.close();verify()
    if sha(database)!=receipt['database_sha256']:raise ValueError('Registry changed during readback')
    result={'status':'passed_full_structure_domain_registry_readback','counts':dict(counts),'candidate_domain_intervals_by_policy':dict(candidates),'producer_receipt_sha256':sha(root/'receipt.json'),'database_sha256':receipt['database_sha256'],'plan_sha256':ph,'elapsed_seconds':time.monotonic()-started,'scope':'Every model, retained interval, policy membership, candidate flag and protein link reconstructed from source data without importing producer interval logic. Overlap and partial-hit exclusions recomputed from intervals and HMM coverage. Shared SQLite sources and SHA helper; does not repeat annotation searches, coordinate validation or establish structural boundaries.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    pins[str(out/'receipt.json')]=sha(out/'receipt.json')
    pins[str(out/'state.json')]=sha(out/'state.json')
    for path,digest in pins.items():assert sha(path)==digest,path
    result.update(source_hashes=pins,scientific_eligibility=False,gpu=False,new_predictions=0,
                  original_producer_transport_sha256=sha(transport_path))
    with args.receipt.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2),flush=True)


if __name__=='__main__':main()
