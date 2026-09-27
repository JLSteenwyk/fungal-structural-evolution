#!/usr/bin/env python3
"""Independently verify every background model join and work partition."""
import argparse,csv,hashlib,json,sqlite3,time
from collections import Counter
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def load_models(path):
    result={}
    with Path(path).open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],m['version']);assert key not in result;result[key]=m
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);a=ap.parse_args();config=json.loads(a.config.read_text());ch=sha(a.config)
    def verify_config():
        assert sha(a.config)==ch
        for p,h in config['pins'].items():assert sha(p)==h,p
    verify_config();dep=config['dependency']
    while True:
        try:
            proc=psutil.Process(dep['pid'])
            if proc.create_time()!=dep['created'] or proc.status()==psutil.STATUS_ZOMBIE:break
            assert proc.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify_config();plan=json.loads(Path(config['plan']).read_text());root=Path(plan['output']);source=Path(plan['source']);r=json.loads((root/'receipt.json').read_text());sr=json.loads((source/'receipt.json').read_text());audit=json.loads(Path(plan['readback']).read_text())
    assert r['status']=='complete_background_measurement_inventory_pending_readback' and r['plan_sha256']==sha(config['plan'])
    assert audit['status']=='passed_full_background_orthology_membership_readback'
    assert r['source_receipt_sha256']==audit['producer_receipt_sha256']==sha(source/'receipt.json') and r['source_readback_sha256']==sha(plan['readback'])
    bindings=dict(plan['pins']);bindings.update({str(root/name):h for name,h in r['artifacts'].items()})
    for path in [root/'receipt.json',source/'receipt.json',Path(plan['readback'])]:bindings[str(path)]=sha(path)
    source_table=source/'candidate_orthology_membership.tsv';bindings[str(source_table)]=sr['artifacts'][source_table.name]
    def verify():
        verify_config()
        for p,h in bindings.items():assert sha(p)==h,p
    verify()
    with source_table.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==r['candidates']==audit['candidate_rows']
    genes={row[k] for row in rows for k in ['gene_a','gene_b']};links={}
    con=sqlite3.connect('file:'+str(Path(plan['bridge']).resolve())+'?mode=ro',uri=True)
    for t,p,seq,m,v,path in con.execute('SELECT taxon_id,protein_id,sequence_sha256,model_id,version,model_path FROM structures'):
        gene=t+'_'+p
        if gene in genes:assert gene not in links;links[gene]=(m,v,seq,path)
    con.close();assert set(links)==genes
    wanted={(v[0],v[1]) for v in links.values()};models={}
    with (Path(plan['catalog'])/'models.jsonl').open() as f:
        for line in f:
            m=json.loads(line);key=(m['model_id'],m['version'])
            if key in wanted:assert key not in models;models[key]=m
    assert set(models)==wanted
    for m,v,seq,path in links.values():assert models[m,v]['sequence_sha256']==seq and models[m,v]['path']==path
    assert load_models(root/'models.jsonl')==models
    existing_pairs={};old_models=set()
    for entry in plan['existing_queues']:
        with (Path(entry['path'])/'model_pairs.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                identity=tuple(sorted([(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]));existing_pairs.setdefault(identity,entry['label'])
        old=load_models(Path(entry['path'])/'models.jsonl');old_models.update(old)
        for key in wanted & old.keys():assert old[key]==models[key]
    counts=Counter();active=set();pairs={};seen=set()
    with (root/'candidate_measurement_links.tsv').open() as f:
        exported=iter(csv.DictReader(f,delimiter='\t'))
        for source_row in rows:
            actual=next(exported,None);assert actual is not None
            identity=(source_row['gene_a'],source_row['gene_b']);assert identity not in seen;seen.add(identity)
            expected=dict(source_row);eligible=[];oriented=[]
            for guide in ['profile','mafft']:
                qualified=source_row[guide+'_native_ortholog']=='1' and source_row[guide+'_candidate_status']=='cross_taxon_unreported_candidate';eligible.append(qualified)
                expected[guide+'_candidate_and_native_ortholog']=str(int(qualified))
            flags=dict(candidate_and_native_ortholog_either=int(sum(eligible)>0),candidate_and_native_ortholog_both=int(sum(eligible)==2),both_guides_and_parents_unreported=int(sum(eligible)==2 and source_row['both_parents_unreported']=='1'))
            expected.update({k:str(v) for k,v in flags.items()});counts.update(flags)
            for suffix in ['a','b']:
                m,v,seq,path=links[source_row['gene_'+suffix]];key=(m,v);oriented.append(key)
                for field in ['model_id','version','sequence_sha256','length','mean_ca_plddt','fraction_ca_plddt_below50']:expected[field+'_'+suffix]=str(models[key][field])
            ordered=tuple(sorted(oriented));pair_key=hashlib.sha256(json.dumps(ordered,separators=(',',':')).encode()).hexdigest();expected['pair_key']=pair_key
            if not any(eligible):disposition='neither_guide_eligible'
            elif len(set(oriented))==1:disposition='identical_model_no_alignment';active.update(oriented)
            else:
                disposition=existing_pairs.get(ordered,'new_model_pair');active.update(oriented)
                pairs[pair_key]=dict(pair_key=pair_key,model_a=ordered[0][0],version_a=str(ordered[0][1]),model_b=ordered[1][0],version_b=str(ordered[1][1]),work_disposition=disposition)
            expected['measurement_disposition']=disposition;counts[disposition]+=1
            assert actual==expected,identity
        assert next(exported,None) is None
    with (root/'model_pairs.tsv').open() as f:
        got=list(csv.DictReader(f,delimiter='\t'));assert len(got)==len(pairs) and {row['pair_key']:row for row in got}==pairs
    assert load_models(root/'active_models.jsonl')=={k:models[k] for k in active}
    additional=active-old_models
    assert load_models(root/'additional_models.jsonl')=={k:models[k] for k in additional}
    assert dict(counts)==r['counts'] and len(models)==r['all_models'] and len(active)==r['active_models'] and len(additional)==r['additional_models'] and len(pairs)==r['distinct_eligible_model_pairs']
    new=sum(row['work_disposition']=='new_model_pair' for row in pairs.values());assert new==r['new_model_pairs']
    verify()
    result=dict(status='passed_full_background_measurement_inventory_readback',config_sha256=ch,producer_receipt_sha256=sha(root/'receipt.json'),candidates=len(rows),all_models=len(models),active_models=len(active),additional_models=len(additional),distinct_eligible_model_pairs=len(pairs),new_model_pairs=new,counts=dict(counts),checker_sha256=sha(__file__),scope='Every candidate field, native eligibility flag, exact bridge/catalog join and metadata value verified; full pair/model output sets and existing/new partitions reconstructed independently. No raw coordinates or alignments validated and no matched controls selected.')
    out=Path(config['output']);assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
