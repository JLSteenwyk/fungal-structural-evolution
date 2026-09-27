#!/usr/bin/env python3
"""End-to-end completed-input fixture: one pair, both orders and two masks."""
import csv,hashlib,json,math,subprocess,sys,tempfile
from pathlib import Path
from duplication_alignment_inputs import render_ca
from run_ortholog_pair_guide_comparison import sha
with tempfile.TemporaryDirectory() as td:
    root=Path(td);queue=root/'queue';inputs=root/'inputs';queue.mkdir();inputs.mkdir()
    def write(p,r):p.write_text(json.dumps(r)+'\n')
    a=('A',1);b=('B',1);pair=hashlib.sha256(json.dumps([a,b],separators=(',',':')).encode()).hexdigest()
    seq='ACDEFGHIKLMNPQRSTVWY';blob,_,_=render_ca(dict(status='validated',sequence=seq,ca_xyz=[[3*math.cos(i),3*math.sin(i),i] for i in range(len(seq))],ca_plddt=[90]*len(seq)))
    rows=[]
    for name in ['A','B']:
        p=inputs/(name+'.pdb');p.write_bytes(blob)
        for mask in ['full','plddt70']:
            row=dict(model_id=name,version=1,mask=mask,status='ready',path=str(p),sha256=sha(p),sequence=seq)
            if name=='B' and mask=='plddt70':row=dict(model_id=name,version=1,mask=mask,status='too_few_retained_residues')
            rows.append(row)
    (inputs/'inputs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows));(queue/'models.jsonl').write_text(''.join(json.dumps({'model_id':n,'version':1})+'\n' for n in ['A','B']))
    with (queue/'model_pairs.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['pair_key','model_a','version_a','model_b','version_b']);w.writerow([pair,'A',1,'B',1])
    write(queue/'receipt.json',{'unique_distinct_model_pairs':1,'artifacts':{n:sha(queue/n) for n in ['models.jsonl','model_pairs.tsv']}})
    ip=root/'input_plan.json';write(ip,{'fixture':True});write(inputs/'receipt.json',{'status':'complete_duplication_alignment_input_materialization','plan_sha256':sha(ip),'queue_receipt_sha256':sha(queue/'receipt.json'),'artifacts':{'inputs.jsonl':sha(inputs/'inputs.jsonl')},'counts':{'full:ready':2,'plddt70:ready':1,'plddt70:too_few_retained_residues':1}})
    plan={'producer':{'pid':99999999,'created':0,'cmdline':[]},'inputs':str(inputs),'input_plan':str(ip),'queue':str(queue),'output':str(root/'out'),'workers':2,'minimum_free_disk_gib':0,'usalign':'/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign','options':['-mol','prot','-mm','0','-outfmt','0','-ter','2'],'per_pair_timeout_seconds':10,'pins':{str(ip):sha(ip)}}
    pp=root/'plan.json';write(pp,plan);subprocess.run([sys.executable,'scripts/run_duplication_alignments.py','--plan',str(pp)],check=True,capture_output=True,text=True)
    out=Path(plan['output']);r=json.loads((out/'receipt.json').read_text());assert r['directed_dispositions']==4 and r['counts']=={'full:aligned':2,'plddt70:input_unavailable':2}
    with (out/'checkpoint_manifest.tsv').open() as f:
        entries=list(csv.DictReader(f,delimiter='\t'));assert len(entries)==4
        for row in entries:assert sha(out/row['path'])==row['sha256']
print('Full runner handoff passed: both native input orders aligned, both short-mask orders retained, all four checkpoint hashes verified.')
