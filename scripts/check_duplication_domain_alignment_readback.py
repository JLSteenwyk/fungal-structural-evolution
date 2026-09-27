#!/usr/bin/env python3
"""Complete native domain alignment and independent numeric readback fixture."""
import csv,hashlib,json,math,subprocess,sys,tempfile
from collections import Counter
from pathlib import Path
from duplication_alignment_inputs import render_ca
from duplication_domain_alignment_handoff import load_handoff
from run_duplication_domain_alignments import run_job
from run_ortholog_pair_guide_comparison import sha

def write(p,r):p.write_text(json.dumps(r)+'\n')
with tempfile.TemporaryDirectory() as td:
    root=Path(td);inventory=root/'inventory';folder=root/'inputs';inventory.mkdir();folder.mkdir()
    intervals=[dict(interval_id=n,model_id='model'+n,version=1,start=1,end=20,original_length=20,source_sha256='raw') for n in ['a','b']]
    (inventory/'intervals.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in intervals))
    pair=hashlib.sha256(json.dumps(['a','b'],separators=(',',':')).encode()).hexdigest()
    (inventory/'domain_pairs.tsv').write_text('domain_pair_key\tinterval_a\tinterval_b\n'+pair+'\ta\tb\n')
    write(inventory/'receipt.json',dict(unique_intervals=2,unique_interval_pairs=1,artifacts={n:sha(inventory/n) for n in ['intervals.jsonl','domain_pairs.tsv']}))
    seq='ACDEFGHIKLMNPQRSTVWY';blob,_,_=render_ca(dict(status='validated',sequence=seq,ca_xyz=[[3*math.cos(i),3*math.sin(i),i] for i in range(20)],ca_plddt=[90]*20));rows=[]
    for interval in intervals:
        path=folder/(interval['interval_id']+'.pdb');path.write_bytes(blob)
        for mask in ['full','plddt70']:
            row=dict(interval,mask=mask,status='ready',path=str(path),sha256=sha(path),sequence=seq,original_positions=list(range(1,21)),retained_residues=20)
            if interval['interval_id']=='b' and mask=='plddt70':row['status']='too_few_retained_residues'
            rows.append(row)
    (folder/'inputs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    ip=root/'ip.json';ap=root/'ap.json';write(ip,{});write(ap,{})
    counts=dict(Counter(r['mask']+':'+r['status'] for r in rows))
    write(folder/'receipt.json',dict(status='complete_duplication_domain_input_materialization_pending_readback',plan_sha256=sha(ip),inventory_receipt_sha256=sha(inventory/'receipt.json'),input_dispositions=4,counts=counts,artifacts={'inputs.jsonl':sha(folder/'inputs.jsonl')}))
    audit=root/'audit.json';write(audit,dict(status='passed_full_duplication_domain_input_readback',plan_sha256=sha(ap),producer_receipt_sha256=sha(folder/'receipt.json'),counts=counts))
    plan=dict(producer=dict(pid=99999999,created=0,cmdline=[]),inventory=str(inventory),inputs=str(folder),input_readback=str(audit),input_plan=str(ip),readback_plan=str(ap),output=str(root/'out'),workers=2,minimum_free_disk_gib=0,usalign='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign',options=['-mol','prot','-mm','0','-outfmt','0','-ter','2'],per_pair_timeout_seconds=10,pins={})
    pp=root/'plan.json';write(pp,plan);subprocess.run([sys.executable,'scripts/run_duplication_domain_alignments.py','--plan',str(pp)],check=True)
    r=json.loads((root/'out/receipt.json').read_text());assert r['directed_dispositions']==4 and r['counts']=={'full:aligned':2,'plddt70:input_unavailable':2}
    inputs,pairs,bindings,mh=load_handoff(plan)
    _,status=run_job((pair,'a','b','full',0),inputs,plan,sha(pp),mh);assert status=='aligned'
    with (root/'out/checkpoint_manifest.tsv').open() as f:
        records=list(csv.DictReader(f,delimiter='\t'))
    assert len(records)==4
    for row in records:assert sha(root/'out'/row['path'])==row['sha256']
    rp=root/'numeric-plan.json';write(rp,dict(source_plan=str(pp),producer=dict(pid=99999999,created=0,cmdline=[]),output=str(root/'numeric'),pins={str(pp):sha(pp)}))
    subprocess.run([sys.executable,'scripts/readback_duplication_domain_alignments.py','--plan',str(rp)],check=True)
    nr=json.loads((root/'numeric/receipt.json').read_text());assert nr['directed_dispositions']==4 and nr['numerically_checked_alignments']==2 and nr['counts']==r['counts']
    # Dropping one excluded disposition is not permitted even with an updated manifest hash.
    cm=root/'out/checkpoint_manifest.tsv';lines=cm.read_text().splitlines();cm.write_text('\n'.join(lines[:-1])+'\n');r['artifacts']['checkpoint_manifest.tsv']=sha(cm);write(root/'out/receipt.json',r)
    ap=json.loads(rp.read_text());ap['output']=str(root/'bad');write(rp,ap)
    failed=subprocess.run([sys.executable,'scripts/readback_duplication_domain_alignments.py','--plan',str(rp)],capture_output=True,text=True)
    assert failed.returncode!=0 and 'Incomplete disposition readback' in failed.stderr
print('Full domain numeric readback passed native both-order alignments and excluded masks; omitted disposition rejected despite updated hash.')
