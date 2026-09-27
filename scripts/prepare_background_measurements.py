#!/usr/bin/env python3
"""Prepare model comparisons for the complete background pool with explicit eligibility."""
import argparse,csv,hashlib,json,shutil,sqlite3,time
from collections import Counter
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def eligibility(row):
    flags={g:row[g+'_candidate_status']=='cross_taxon_unreported_candidate' and row[g+'_native_ortholog']=='1' for g in ['profile','mafft']}
    return dict(profile_candidate_and_native_ortholog=int(flags['profile']),mafft_candidate_and_native_ortholog=int(flags['mafft']),
                candidate_and_native_ortholog_either=int(any(flags.values())),
                candidate_and_native_ortholog_both=int(all(flags.values())),
                both_guides_and_parents_unreported=int(all(flags.values()) and row['both_parents_unreported']=='1'))


def pair_identity(left,right):
    models=sorted([left,right]);key=hashlib.sha256(json.dumps(models,separators=(',',':')).encode()).hexdigest()
    return key,models


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();dep=plan['dependency']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();source=Path(plan['source']);sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text());ah=sha(plan['readback']);rh=sha(source/'receipt.json')
    assert sr['status']=='complete_background_native_ortholog_membership_pending_readback'
    assert audit['status']=='passed_full_background_orthology_membership_readback' and audit['producer_receipt_sha256']==rh
    table=source/'candidate_orthology_membership.tsv';assert sha(table)==sr['artifacts'][table.name]
    with table.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==audit['candidate_rows']==sr['candidate_rows']
    needed={r[k] for r in rows for k in ['gene_a','gene_b']};links={}
    db=sqlite3.connect('file:'+str(Path(plan['bridge']).resolve())+'?mode=ro',uri=True)
    for taxon,protein,seq,model,version,path in db.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id,version,model_path FROM structures'):
        gene=taxon+'_'+protein
        if gene in needed:
            assert gene not in links;links[gene]=(seq,model,version,path)
    db.close();assert set(links)==needed
    modelkeys={(v[1],v[2]) for v in links.values()};models={};catalog=Path(plan['catalog']);cr=json.loads((catalog/'receipt.json').read_text())
    assert sha(catalog/'models.jsonl')==cr['artifacts']['models.jsonl']
    with (catalog/'models.jsonl').open() as f:
        for line in f:
            model=json.loads(line);key=(model['model_id'],model['version'])
            if key in modelkeys:assert key not in models;models[key]=model
    assert set(models)==modelkeys
    for seq,m,v,path in links.values():assert (models[m,v]['sequence_sha256'],models[m,v]['path'])==(seq,path)
    previous_pairs={};previous_models=set()
    for entry in plan['existing_queues']:
        root=Path(entry['path']);r=json.loads((root/'receipt.json').read_text())
        for filename in ['model_pairs.tsv','models.jsonl']:assert sha(root/filename)==r['artifacts'][filename]
        with (root/'model_pairs.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                values=[(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))];key,_=pair_identity(*values);assert key==row['pair_key']
                previous_pairs.setdefault(key,entry['label'])
        with (root/'models.jsonl').open() as f:
            for line in f:
                m=json.loads(line);key=(m['model_id'],m['version']);previous_models.add(key)
                if key in models:assert m==models[key]
    out=Path(plan['output']);assert shutil.disk_usage(out.parent).free>=plan['resources']['minimum_free_disk_bytes'];out.mkdir(exist_ok=False)
    counts=Counter();queued={};active=set();writer=None
    with (out/'candidate_measurement_links.tsv').open('w') as f:
        for row in rows:
            record=dict(row,**eligibility(row));oriented=[]
            for suffix in ['a','b']:
                seq,m,v,path=links[row['gene_'+suffix]];oriented.append((m,v));model=models[m,v]
                for guide in ['profile','mafft']:
                    if row[guide+'_model_'+suffix]:assert row[guide+'_model_'+suffix]==m and int(row[guide+'_version_'+suffix])==v
                for field in ['model_id','version','sequence_sha256','length','mean_ca_plddt','fraction_ca_plddt_below50']:
                    record[field+'_'+suffix]=model[field]
            key,ordered=pair_identity(*oriented);record['pair_key']=key
            if not record['candidate_and_native_ortholog_either']:status='neither_guide_eligible'
            elif oriented[0]==oriented[1]:status='identical_model_no_alignment';active.update(oriented)
            else:
                status=previous_pairs.get(key,'new_model_pair');queued[key]=ordered;active.update(oriented)
            record['measurement_disposition']=status;counts[status]+=1
            for flag in ['candidate_and_native_ortholog_either','candidate_and_native_ortholog_both','both_guides_and_parents_unreported']:counts[flag]+=record[flag]
            if writer is None:writer=csv.DictWriter(f,list(record),delimiter='\t',lineterminator='\n');writer.writeheader()
            writer.writerow(record)
    with (out/'model_pairs.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['pair_key','model_a','version_a','model_b','version_b','work_disposition'])
        for key,(x,y) in sorted(queued.items()):w.writerow([key,*x,*y,previous_pairs.get(key,'new_model_pair')])
    for filename,keys in [('models.jsonl',modelkeys),('active_models.jsonl',active),('additional_models.jsonl',active-previous_models)]:
        with (out/filename).open('w') as f:
            for key in sorted(keys):f.write(json.dumps(models[key],separators=(',',':'))+'\n')
    verify();assert sha(plan['readback'])==ah and sha(source/'receipt.json')==rh and sha(table)==sr['artifacts'][table.name]
    result=dict(status='complete_background_measurement_inventory_pending_readback',plan_sha256=ph,source_receipt_sha256=rh,source_readback_sha256=ah,candidates=len(rows),counts=dict(counts),all_models=len(models),active_models=len(active),additional_models=len(active-previous_models),distinct_eligible_model_pairs=len(queued),new_model_pairs=sum(k not in previous_pairs for k in queued),artifacts={p.name:sha(p) for p in out.iterdir()},scope='All modeled background candidates retained with native orthology and guide/parent flags, exact frozen model joins, length/confidence covariates and work partitions. Existing queue membership does not prove completed alignment or valid coordinates. Identical models remain explicit; neither-guide eligibility retained. Not matched controls; domain/coverage/phylogenetic matching and coordinate/readback checks remain required.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
