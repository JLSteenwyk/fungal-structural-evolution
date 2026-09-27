#!/usr/bin/env python3
"""Exercise complete background handoff with masks, rejected sources and corrupt proof."""
import gzip,json,subprocess,sys,tempfile
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha
with tempfile.TemporaryDirectory() as td:
    root=Path(td);coords=root/'coords';audit=root/'audit';inv=root/'inventory'
    for p in [coords,audit,inv]:p.mkdir()
    def write(p,x):p.write_text(json.dumps(x)+'\n')
    rows=[]
    for ident,conf in [('A',[70,20,90,80,60]),('B',[0]*5)]:rows.append(dict(model_id=ident,version=1,status='validated',sequence='ACDEF',length=5,ca_xyz=[[i,i/2,0] for i in range(5)],ca_plddt=conf,source_sha256='fixture'))
    rows.append(dict(model_id='C',version=1,status='rejected_content',reason='fixture rejection',source_sha256='fixture'))
    for row in rows:row.update(source_path='fixture.cif',sequence_sha256='fixture-sequence')
    path=coords/'shard-00000.jsonl.gz'
    with gzip.open(path,'wt') as f:
        for row in rows:f.write(json.dumps(row)+'\n')
    models=[dict(model_id=r['model_id'],version=1,path=r['source_path'],sha256=r['source_sha256'],sequence_sha256=r['sequence_sha256']) for r in rows]
    (inv/'additional_models.jsonl').write_text(''.join(json.dumps(m)+'\n' for m in models))
    write(inv/'receipt.json',dict(status='complete_background_measurement_inventory_pending_readback',additional_models=3,artifacts={'additional_models.jsonl':sha(inv/'additional_models.jsonl')}))
    ia=root/'inventory-audit.json';write(ia,dict(status='passed_full_background_measurement_inventory_readback',producer_receipt_sha256=sha(inv/'receipt.json'),additional_models=3))
    cp=root/'coordinate-plan.json';rp=root/'readback-plan.json'
    write(cp,dict(inventory=str(inv),readback=str(ia),output=str(coords)))
    write(rp,dict(source_plan=str(cp),source=str(coords),models=str(inv/'additional_models.jsonl')))
    write(coords/'receipt.json',dict(status='complete_background_coordinate_validation_with_dispositions',models=3,plan_sha256=sha(cp),source_inventory_receipt_sha256=sha(inv/'receipt.json'),shards=[dict(output=path.name,output_sha256=sha(path))]))
    proof=audit/(path.name+'.readback.json');write(proof,dict(source_sha256=sha(path)))
    write(audit/'receipt.json',dict(status='passed_background_exported_ca_readback',plan_sha256=sha(rp),producer_receipt_sha256=sha(coords/'receipt.json'),models=3,proofs={proof.name:sha(proof)}))
    plan=dict(producer=dict(pid=99999999,created=0,cmdline=[]),readback=str(audit),readback_plan=str(rp),coordinates=str(coords),coordinate_plan=str(cp),output=str(root/'out'),minimum_free_disk_gib=0,pins={str(cp):sha(cp),str(rp):sha(rp)})
    pp=root/'plan.json';write(pp,plan)
    subprocess.run([sys.executable,'scripts/materialize_background_inputs.py','--plan',str(pp)],check=True,capture_output=True,text=True)
    out=Path(plan['output']);receipt=json.loads((out/'receipt.json').read_text())
    assert receipt['status']=='complete_background_alignment_input_materialization' and receipt['models']==3 and receipt['input_dispositions']==6
    assert receipt['counts']=={'full:ready':2,'plddt70:ready':1,'plddt70:too_few_retained_residues':1,'full:source_rejected':1,'plddt70:source_rejected':1}
    entries=[json.loads(l) for l in (out/'inputs.jsonl').read_text().splitlines()]
    selected=next(x for x in entries if x['model_id']=='A' and x['mask']=='plddt70');assert selected['sequence']=='ADE' and selected['original_positions']==[1,3,4]
    for entry in entries:
        if entry['status']=='ready':assert sha(entry['path'])==entry['sha256']
    proof.write_text('{}');plan['output']=str(root/'bad');write(pp,plan)
    bad=subprocess.run([sys.executable,'scripts/materialize_background_inputs.py','--plan',str(pp)],capture_output=True,text=True)
    assert bad.returncode!=0 and 'Unaudited or changed coordinate shard' in bad.stderr
print('Complete handoff, six dispositions, sparse original positions, written hashes and changed-proof rejection passed.')
