#!/usr/bin/env python3
"""Synthetic completed-handoff fixture for full, masked, short and rejected inputs."""
import gzip,json,subprocess,sys,tempfile
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha

with tempfile.TemporaryDirectory() as td:
    root=Path(td);coords=root/'coords';audit=root/'audit';queue=root/'queue'
    for p in [coords,audit,queue]:p.mkdir()
    def write(path,obj):path.write_text(json.dumps(obj)+'\n')
    rows=[]
    for ident,confidence in [('A',[70,20,90,80,60]),('B',[0]*5)]:rows.append(dict(model_id=ident,version=1,status='validated',sequence='ACDEF',length=5,ca_xyz=[[i,i/2,0] for i in range(5)],ca_plddt=confidence,source_sha256='fixture'))
    rows.append(dict(model_id='C',version=1,status='rejected_content',reason='fixture rejection',source_sha256='fixture'))
    path=coords/'shard-00000.jsonl.gz'
    with gzip.open(path,'wt') as f:
        for row in rows:f.write(json.dumps(row)+'\n')
    coordinate_plan=root/'coordinate_plan.json';readback_plan=root/'readback_plan.json';write(coordinate_plan,{'queue':str(queue)});write(readback_plan,{'fixture':True})
    write(queue/'receipt.json',{'status':'complete_reviewed_duplication_model_pair_queue','unique_models':3})
    source_receipt={'status':'complete_duplication_coordinate_validation_with_dispositions','models':3,'plan_sha256':sha(coordinate_plan),'shards':[{'output':path.name,'output_sha256':sha(path)}]};write(coords/'receipt.json',source_receipt)
    proof=audit/(path.name+'.readback.json');write(proof,{'source_sha256':sha(path)})
    write(audit/'receipt.json',{'status':'passed_duplication_exported_ca_readback','plan_sha256':sha(readback_plan),'producer_receipt_sha256':sha(coords/'receipt.json'),'models':3,'proofs':{proof.name:sha(proof)}})
    plan={'producer':{'pid':99999999,'created':0,'cmdline':[]},'readback':str(audit),'readback_plan':str(readback_plan),'coordinates':str(coords),'coordinate_plan':str(coordinate_plan),'queue':str(queue),'output':str(root/'output'),'minimum_free_disk_gib':0,'pins':{str(coordinate_plan):sha(coordinate_plan),str(readback_plan):sha(readback_plan)}}
    pp=root/'plan.json';write(pp,plan)
    subprocess.run([sys.executable,'scripts/materialize_duplication_alignment_inputs.py','--plan',str(pp)],check=True,capture_output=True,text=True)
    out=Path(plan['output']);r=json.loads((out/'receipt.json').read_text());assert r['models']==3 and r['input_dispositions']==6
    assert r['counts']=={'full:ready':2,'plddt70:ready':1,'plddt70:too_few_retained_residues':1,'full:source_rejected':1,'plddt70:source_rejected':1}
    entries=[json.loads(l) for l in (out/'inputs.jsonl').read_text().splitlines()]
    m=next(x for x in entries if x['model_id']=='A' and x['mask']=='plddt70');assert m['sequence']=='ADE' and m['original_positions']==[1,3,4]
    for e in entries:
        if e['status']=='ready':assert sha(e['path'])==e['sha256']
    # Corrupted audit proof must stop a fresh output before coordinate use.
    write(proof,{'source_sha256':'wrong'});plan['output']=str(root/'bad_output');write(pp,plan)
    result=subprocess.run([sys.executable,'scripts/materialize_duplication_alignment_inputs.py','--plan',str(pp)],capture_output=True,text=True);assert result.returncode!=0 and 'Unaudited or changed coordinate shard' in result.stderr
print('Completed-handoff fixture passed all six dispositions, written hashes and sparse residue mapping; altered audit proof rejected.')
