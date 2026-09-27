#!/usr/bin/env python3
"""Exercise three audited partitions, pair reuse, native masks/orders and corruption gates."""
import csv,hashlib,json,math,subprocess,sys,tempfile
from collections import Counter
from pathlib import Path
from duplication_alignment_inputs import render_ca
from background_alignment_handoff import load_handoff
from run_ortholog_pair_guide_comparison import sha

def write(p,x):p.write_text(json.dumps(x)+'\n')
def lines(p,rows):p.write_text(''.join(json.dumps(x)+'\n' for x in rows))
def pair(a,b,d):
    return dict(pair_key=hashlib.sha256(json.dumps([(a,1),(b,1)],separators=(',',':')).encode()).hexdigest(),model_a=a,version_a=1,model_b=b,version_b=1,work_disposition=d)
def table(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with tempfile.TemporaryDirectory() as td:
    root=Path(td);inv=root/'inventory';inv.mkdir();models=[dict(model_id=n,version=1,sha256='source-'+n) for n in 'ABC']
    lines(inv/'active_models.jsonl',models)
    table(inv/'model_pairs.tsv',[pair('A','B','already_in_reference_queue'),pair('A','C','new_model_pair'),pair('B','C','new_model_pair')])
    write(inv/'receipt.json',dict(status='complete_background_measurement_inventory_pending_readback',distinct_eligible_model_pairs=3,new_model_pairs=2,artifacts={n:sha(inv/n) for n in ['active_models.jsonl','model_pairs.tsv']}))
    proof=root/'proof.json';write(proof,dict(status='passed_full_background_measurement_inventory_readback',producer_receipt_sha256=sha(inv/'receipt.json'),active_models=3))
    seq='ACDEFGHIKLMNPQRSTVWY';blob,_,_=render_ca(dict(status='validated',sequence=seq,ca_xyz=[[3*math.cos(i),3*math.sin(i),i] for i in range(len(seq))],ca_plddt=[90]*len(seq)))
    sources={};queues=[]
    for m,label in zip(models,['primary','reference','background']):
        source=root/(label+'-models');source.mkdir();lines(source/'models.jsonl',[m])
        table(source/'model_pairs.tsv',[pair('A','B','already_in_reference_queue')])
        write(source/'receipt.json',dict(artifacts={n:sha(source/n) for n in ['models.jsonl','model_pairs.tsv']}))
        if label=='reference':queues.append(dict(path=str(source),label='already_in_reference_queue'))
        folder=root/label;folder.mkdir();pdb=folder/'input.pdb';pdb.write_bytes(blob)
        coord=root/(label+'-coords');coord.mkdir();write(coord/'receipt.json',dict(status='validated'))
        audit=root/(label+'-audit');audit.mkdir();write(audit/'receipt.json',dict(status='audited',producer_receipt_sha256=sha(coord/'receipt.json')))
        ip=root/(label+'-plan.json');write(ip,dict(output=str(folder),readback=str(audit),coordinates=str(coord)))
        rows=[dict(model_id=m['model_id'],version=1,source_sha256=m['sha256'],mask=mask,status=('too_few_retained_residues' if label=='background' and mask=='plddt70' else 'ready'),path=str(pdb),sha256=sha(pdb),sequence=seq) for mask in ['full','plddt70']]
        lines(folder/'inputs.jsonl',rows)
        write(folder/'receipt.json',dict(status='materialized',plan_sha256=sha(ip),inventory_receipt_sha256=sha(source/'receipt.json'),readback_receipt_sha256=sha(audit/'receipt.json'),coordinate_receipt_sha256=sha(coord/'receipt.json'),models=1,input_dispositions=2,counts=dict(Counter(r['mask']+':'+r['status'] for r in rows)),artifacts={'inputs.jsonl':sha(folder/'inputs.jsonl')}))
        sources[label]=dict(inputs=str(folder),input_plan=str(ip),model_inventory=str(source),model_file='models.jsonl',status='materialized',inventory_receipt_field='inventory_receipt_sha256',coordinate_readback_status='audited')
    plan=dict(inventory=str(inv),inventory_readback=str(proof),input_sources=sources,existing_queues=queues,resources=dict(new_model_pairs=2),output=str(root/'out'),workers=2,minimum_free_disk_gib=0,usalign='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign',options=['-mol','prot','-mm','0','-outfmt','0','-ter','2'],per_pair_timeout_seconds=10,pins={})
    inputs,pairs,_,_=load_handoff(plan);assert len(inputs)==6 and len(pairs)==2
    pp=root/'plan.json';write(pp,plan);subprocess.run([sys.executable,'scripts/run_background_alignments.py','--plan',str(pp)],check=True)
    r=json.loads((root/'out/receipt.json').read_text());assert r['directed_dispositions']==8 and r['counts']=={'full:aligned':4,'plddt70:input_unavailable':4}
    with (root/'out/checkpoint_manifest.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):assert sha(root/'out'/r['path'])==r['sha256']
    def reject(phrase):
        try:load_handoff(plan)
        except ValueError as e:assert phrase in str(e),str(e)
        else:raise AssertionError('Corruption accepted: '+phrase)
    ap=root/'background-audit/receipt.json';original=ap.read_text();write(ap,dict(status='wrong'));reject('Coordinate proof');ap.write_text(original)
    manifest=root/'background/inputs.jsonl';original=manifest.read_text();rp=root/'background/receipt.json';rr=json.loads(rp.read_text())
    rows=[json.loads(l) for l in original.splitlines()];rows[0]['source_sha256']='wrong';lines(manifest,rows);rr['artifacts']['inputs.jsonl']=sha(manifest);write(rp,rr);reject('Wrong input source')
    rows=[json.loads(l) for l in original.splitlines()];lines(manifest,rows[:1]);rr['artifacts']['inputs.jsonl']=sha(manifest);write(rp,rr);reject('Incomplete input grid')
print('Passed three-source native handoff, existing-pair exclusion, eight dispositions, checkpoint hashes, changed-proof/source and missing-mask rejection.')
