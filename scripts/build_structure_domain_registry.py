#!/usr/bin/env python3
"""Join the full structural catalog to policy-preserving annotation intervals."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import json
from pathlib import Path
import shutil
import sqlite3
import time
from catalog_whole_proteome_structures import sha

POLICIES={'alignment_evalue','alignment_bitscore','envelope_evalue','envelope_bitscore'}


def intervals(payload,length,agree):
    if set(payload['policies'])!=POLICIES:raise ValueError('Unexpected policy grid')
    hits={};members=[]
    for policy,entry in sorted(payload['policies'].items()):
        alt=payload['alternatives'][entry['alternative_index']]
        conservative=bool(agree) and not any([entry['primary_rank_tie_pairs'],entry['unresolved_overlap_pairs'],entry['candidate_nested_pairs'],alt['alignment_overlap_pairs'],alt['partial_hmm_hits_below_070']])
        seen=set()
        for hit in alt['annotations']:
            hid=hit['hit_id']
            if hid in seen:raise ValueError('Repeated hit in policy')
            seen.add(hid)
            start,end,es,ee=[hit[k] for k in ('alignment_start','alignment_end','envelope_start','envelope_end')]
            if not all(isinstance(x,int) for x in (start,end,es,ee)) or not 1<=es<=start<=end<=ee<=length:raise ValueError('Invalid sequence interval')
            coverage=Decimal(hit['hmm_coverage'])
            if not coverage.is_finite() or not 0<=coverage<=1:raise ValueError('Invalid HMM coverage')
            row=(hid,hit['pfam_accession'],hit['pfam_type'],hit['pfam_clan'],start,end,es,ee,str(coverage))
            if hid in hits and hits[hid]!=row:raise ValueError('Inconsistent hit across policies')
            hits[hid]=row
            eligible=conservative and hit['pfam_type']=='Domain' and coverage>=Decimal('0.70') and end-start+1>=30
            members.append((policy,hid,int(conservative),int(eligible)))
    return list(hits.values()),members


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Pinned source changed: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient disk')
    out.mkdir(parents=True);dbpath=out/'structure_domains.sqlite'
    db=sqlite3.connect(dbpath);db.execute('PRAGMA foreign_keys=ON')
    db.executescript('''CREATE TABLE models(model_key TEXT PRIMARY KEY,model_id TEXT,version INTEGER,sequence_sha256 TEXT UNIQUE,length INTEGER,path TEXT,raw_hits INTEGER,policies_agree INTEGER);
CREATE TABLE segments(model_key TEXT,hit_id TEXT,pfam_accession TEXT,pfam_type TEXT,pfam_clan TEXT,alignment_start INTEGER,alignment_end INTEGER,envelope_start INTEGER,envelope_end INTEGER,hmm_coverage TEXT,PRIMARY KEY(model_key,hit_id),FOREIGN KEY(model_key) REFERENCES models(model_key));
CREATE TABLE policy_hits(model_key TEXT,policy TEXT,hit_id TEXT,conservative_architecture INTEGER,domain_interval_candidate INTEGER,PRIMARY KEY(model_key,policy,hit_id),FOREIGN KEY(model_key,hit_id) REFERENCES segments(model_key,hit_id));
CREATE TABLE protein_links(taxon_id TEXT,protein_id TEXT,model_key TEXT,PRIMARY KEY(taxon_id,protein_id),FOREIGN KEY(model_key) REFERENCES models(model_key));''')
    annotations=sqlite3.connect('file:'+str(Path(plan['architectures']).resolve())+'?mode=ro',uri=True)
    counts=Counter();policy_counts=Counter();started=time.monotonic()
    with Path(plan['models']).open() as handle:
        for line in handle:
            m=json.loads(line);key=Path(m['path']).stem
            found=annotations.execute('SELECT raw_hits,retained_sets_agree,candidate_architectures_json FROM queries WHERE sequence_id=?',('S'+m['sequence_sha256'],)).fetchone()
            if found is None:raise ValueError('Model sequence absent from full annotation inventory')
            raw,agree,text=found;hits,members=intervals(json.loads(text),m['length'],agree)
            if bool(raw)!=bool(hits):raise ValueError('Raw-hit/retained-hit disposition inconsistent')
            db.execute('INSERT INTO models VALUES (?,?,?,?,?,?,?,?)',(key,m['model_id'],m['version'],m['sequence_sha256'],m['length'],m['path'],raw,agree))
            db.executemany('INSERT INTO segments VALUES (?,?,?,?,?,?,?,?,?,?)',[(key,*r) for r in hits])
            db.executemany('INSERT INTO policy_hits VALUES (?,?,?,?,?)',[(key,*r) for r in members])
            counts['models']+=1;counts['models_without_hits']+=not raw;counts['segments']+=len(hits);counts['policy_hits']+=len(members)
            for policy,_,_,eligible in members:policy_counts[policy]+=eligible
            if counts['models']%10000==0:
                db.commit();(out/'state.json').write_text(json.dumps({'stage':'joining_models','counts':dict(counts),'elapsed_seconds':time.monotonic()-started})+'\n')
    with Path(plan['protein_links']).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=Path(row['model_path']).stem
            model=db.execute('SELECT sequence_sha256,model_id,version FROM models WHERE model_key=?',(key,)).fetchone()
            if model!=(row['sequence_sha256'],row['model_id'],int(row['version'])):raise ValueError('Catalog protein/model identity differs')
            annotation=annotations.execute('SELECT sequence_id FROM proteins WHERE taxon_id=? AND protein_id=?',(row['taxon_id'],row['protein_id'])).fetchone()
            if annotation!=('S'+row['sequence_sha256'],):raise ValueError('Annotation protein/sequence identity differs')
            db.execute('INSERT INTO protein_links VALUES (?,?,?)',(row['taxon_id'],row['protein_id'],key));counts['protein_links']+=1
    db.commit()
    if counts['models']!=plan['expected_models'] or counts['protein_links']!=plan['expected_protein_links']:raise ValueError('Catalog scope differs')
    if db.execute('PRAGMA integrity_check').fetchone()!=('ok',) or db.execute('PRAGMA foreign_key_check').fetchall():raise ValueError('Database integrity failure')
    for table,count in [('models','models'),('segments','segments'),('policy_hits','policy_hits'),('protein_links','protein_links')]:
        if db.execute('SELECT count(*) FROM '+table).fetchone()[0]!=counts[count]:raise ValueError('Database count mismatch')
    db.close();annotations.close();verify()
    receipt={'status':'complete_structure_domain_interval_registry_pending_independent_readback','counts':dict(counts),'candidate_domain_intervals_by_policy':dict(policy_counts),'plan_sha256':ph,'database_sha256':sha(dbpath),'elapsed_seconds':time.monotonic()-started,'coordinate_convention':'1-based inclusive source protein positions; alignment and envelope boundaries retained separately','scope':'Full exact-sequence model/annotation and protein-link join with all policy alternatives and Pfam types. Candidate Domain intervals require conservative architecture, HMM coverage>=0.70 and alignment length>=30. No residue confidence or PAE qualification, coordinate extraction, validated structural boundaries, homology or evolutionary events.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
