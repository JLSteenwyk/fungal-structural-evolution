#!/usr/bin/env python3
"""End-to-end domain serialization readback and numeric corruption rejection fixture."""
import gzip,io,json,subprocess,sys,tempfile
from pathlib import Path
from Bio.PDB import PDBParser
from duplication_domain_inputs import render_interval
from run_ortholog_pair_guide_comparison import sha

record=dict(status='validated',sequence='ACDEFGHIKL',length=10,ca_xyz=[[i+.12345,i/2,0] for i in range(10)],ca_plddt=[90,90,90,20,90,90,20,90,90,90])
for threshold,expected in [(None,[3,4,5,6,7,8]),(70,[3,5,6,8])]:
    blob,seq,positions=render_interval(record,3,8,threshold)
    residues=list(PDBParser(QUIET=True).get_structure('fixture',io.StringIO(blob.decode())).get_residues())
    assert positions==expected==[r.id[1] for r in residues]
    for residue,pos in zip(residues,expected):
        assert max(abs(float(x)-y) for x,y in zip(residue['CA'].coord,record['ca_xyz'][pos-1]))<.000501
        assert residue['CA'].bfactor==record['ca_plddt'][pos-1]
assert render_interval(dict(record,ca_plddt=[0]*10),3,8,70)[2]==[]
for start,end in [(0,3),(4,3),(1,11)]:
    try:render_interval(record,start,end)
    except ValueError:pass
    else:raise AssertionError('Invalid interval accepted')

def write(p,r):p.write_text(json.dumps(r)+'\n')
with tempfile.TemporaryDirectory() as td:
    root=Path(td);inventory=root/'inventory';inventory.mkdir();intervals=[];sources=[]
    for index in [0,1]:
        coords=root/f'coords{index}';audit=root/f'audit{index}';coords.mkdir();audit.mkdir()
        r=dict(record,model_id=str(index),version=1,source_path='fixture',source_sha256='raw',sequence_sha256='seq')
        if index==1:r.update(status='rejected_content',reason='fixture rejection')
        path=coords/'shard-00000.jsonl.gz'
        with gzip.open(path,'wt') as f:f.write(json.dumps(r)+'\n')
        cp=root/f'cp{index}.json';ap=root/f'ap{index}.json';write(cp,dict(output=str(coords)));write(ap,dict(source=str(coords),source_plan=str(cp)))
        write(coords/'receipt.json',dict(status='complete_duplication_coordinate_validation_with_dispositions',plan_sha256=sha(cp),shards=[dict(output=path.name,output_sha256=sha(path))]))
        proof=audit/(path.name+'.readback.json');write(proof,dict(source_sha256=sha(path)))
        write(audit/'receipt.json',dict(status='passed_duplication_exported_ca_readback',plan_sha256=sha(ap),producer_receipt_sha256=sha(coords/'receipt.json'),proofs={proof.name:sha(proof)}))
        sources.append(dict(producer=dict(pid=99999999,created=0,cmdline=[]),readback=str(audit),coordinates=str(coords),readback_plan=str(ap),coordinate_plan=str(cp)))
        intervals.append(dict(interval_id='interval'+str(index),model_id=str(index),version=1,source_path='fixture',source_sha256='raw',sequence_sha256='seq',original_length=10,start=3,end=8,length=6))
    (inventory/'intervals.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in intervals))
    write(inventory/'receipt.json',dict(unique_intervals=2,artifacts={'intervals.jsonl':sha(inventory/'intervals.jsonl')}))
    ir=root/'ir.json';write(ir,dict(status='passed_full_duplication_domain_pair_inventory_readback',producer_receipt_sha256=sha(inventory/'receipt.json')))
    plan=dict(coordinate_sources=sources,inventory=str(inventory),inventory_readback=str(ir),output=str(root/'out'),minimum_free_disk_gib=0,pins={})
    pp=root/'plan.json';write(pp,plan);command=[sys.executable,'scripts/materialize_duplication_domain_inputs.py','--plan',str(pp)]
    subprocess.run(command,check=True,capture_output=True,text=True)
    out=root/'out';r=json.loads((out/'receipt.json').read_text());assert r['input_dispositions']==4 and r['counts']=={'full:ready':1,'plddt70:ready':1,'full:source_rejected':1,'plddt70:source_rejected':1}
    rows=[json.loads(l) for l in (out/'inputs.jsonl').read_text().splitlines()];assert rows[1]['original_positions']==[3,5,6,8]
    for row in rows:
        if row['status']=='ready':assert sha(row['path'])==row['sha256']
    ap=root/'readback-plan.json';write(ap,dict(source_plan=str(pp),producer=dict(pid=99999999,created=0,cmdline=[]),output=str(root/'readback.json'),pins={str(pp):sha(pp)}))
    subprocess.run([sys.executable,'scripts/readback_duplication_domain_inputs.py','--plan',str(ap)],check=True,capture_output=True,text=True)
    audit=json.loads((root/'readback.json').read_text());assert audit['input_dispositions']==4 and audit['verified_pdb_ca_atoms']==10
    # Corrupt a coordinate and update its manifest hash: numeric reconstruction still fails.
    ready=rows[0];pdb=Path(ready['path']);lines=pdb.read_text().splitlines();lines[0]=lines[0][:30]+f'{999.:8.3f}'+lines[0][38:];pdb.write_text('\n'.join(lines)+'\n');ready['sha256']=sha(pdb)
    from readback_duplication_domain_inputs import check_input
    try:check_input(ready,intervals[0],record)
    except ValueError as e:assert 'coordinate differs' in str(e)
    else:raise AssertionError('Altered coordinate accepted')
print('Full two-source domain serialization readback passed; changed coordinate rejected despite updated PDB hash.')
