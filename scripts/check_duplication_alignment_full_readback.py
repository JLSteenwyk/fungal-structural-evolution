#!/usr/bin/env python3
"""End-to-end reference alignment and independent full numeric readback fixture."""
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from duplication_alignment_inputs import render_ca
from duplication_reference_alignment_handoff import load_handoff
from run_ortholog_pair_guide_comparison import sha


def write(path, row):
    path.write_text(json.dumps(row)+'\n')


def pair(a,b,disposition):
    return dict(pair_key=hashlib.sha256(json.dumps([(a,1),(b,1)],separators=(',',':')).encode()).hexdigest(),
                model_a=a,version_a=1,model_b=b,version_b=1,work_disposition=disposition)


def table(path, rows):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


with tempfile.TemporaryDirectory() as directory:
    root=Path(directory);base=root/'base';inventory=root/'inventory'
    base.mkdir();inventory.mkdir()
    models={n:dict(model_id=n,version=1,sha256='source-'+n) for n in 'ABC'}
    for folder,name,names in [(base,'models.jsonl','AB'),(inventory,'models.jsonl','ABC'),
                              (inventory,'additional_models.jsonl','C')]:
        (folder/name).write_text(''.join(json.dumps(models[n])+'\n' for n in names))
    table(base/'model_pairs.tsv',[pair('A','B','existing_duplicate_pair')])
    table(inventory/'model_pairs.tsv',[pair('A','B','existing_duplicate_pair'),pair('A','C','additional_pair')])
    write(base/'receipt.json',dict(status='complete_reviewed_duplication_model_pair_queue',unique_models=2,
          artifacts={name:sha(base/name) for name in ['models.jsonl','model_pairs.tsv']}))
    write(inventory/'receipt.json',dict(status='complete_provisional_reference_comparison_inventory',
          additional_models=1,all_reference_comparison_models=3,additional_model_pairs=1,
          existing_duplicate_model_pairs=1,unique_distinct_model_pairs=2,
          source_pins={str(base/name):sha(base/name) for name in ['receipt.json','models.jsonl']},
          artifacts={name:sha(inventory/name) for name in ['models.jsonl','additional_models.jsonl','model_pairs.tsv']}))
    ledger=root/'ledger.json'
    write(ledger,dict(status='passed_full_reference_comparison_ledger_readback',additional_models=1,
                     unique_models=3,producer_receipt_sha256=sha(inventory/'receipt.json')))
    seq='ACDEFGHIKLMNPQRSTVWY'
    blob,_,_=render_ca(dict(status='validated',sequence=seq,
             ca_xyz=[[3*math.cos(i),3*math.sin(i),i] for i in range(len(seq))],ca_plddt=[90]*len(seq)))
    sources={}
    for label,names in [('primary','AB'),('additional','C')]:
        folder=root/label;folder.mkdir();rows=[]
        for n in names:
            pdb=folder/(n+'.pdb');pdb.write_bytes(blob)
            for mask in ['full','plddt70']:
                row=dict(model_id=n,version=1,mask=mask,source_sha256=models[n]['sha256'],
                         status='ready',path=str(pdb),sha256=sha(pdb),sequence=seq,original_positions=list(range(1,len(seq)+1)),retained_residues=len(seq))
                if n=='C' and mask=='plddt70':row['status']='too_few_retained_residues'
                rows.append(row)
        (folder/'inputs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
        ip=root/(label+'-plan.json');write(ip,dict(output=str(folder)))
        r=dict(status=('complete_duplication_alignment_input_materialization' if label=='primary' else
                       'complete_additional_reference_alignment_input_materialization'),
               plan_sha256=sha(ip),models=len(names),input_dispositions=len(rows),
               counts=dict(Counter(r['mask']+':'+r['status'] for r in rows)),
               artifacts={'inputs.jsonl':sha(folder/'inputs.jsonl')})
        r['queue_receipt_sha256' if label=='primary' else 'inventory_receipt_sha256']=sha((base if label=='primary' else inventory)/'receipt.json')
        write(folder/'receipt.json',r)
        sources[label]=dict(inputs=str(folder),input_plan=str(ip),producer=dict(pid=99999999,created=0,cmdline=[]))
    plan=dict(inventory=str(inventory),base_queue=str(base),inventory_readback=str(ledger),input_sources=sources,
              output=str(root/'out'),workers=2,minimum_free_disk_gib=0,
              usalign='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign',
              options=['-mol','prot','-mm','0','-outfmt','0','-ter','2'],per_pair_timeout_seconds=10,pins={})
    inputs,pairs,bindings,bundle=load_handoff(plan)
    assert len(inputs)==6 and len(pairs)==1 and pairs[0]['model_b']=='C'
    pp=root/'plan.json';write(pp,plan)
    subprocess.run([sys.executable,'scripts/run_duplication_reference_alignments.py','--plan',str(pp)],check=True)
    out=Path(plan['output']);r=json.loads((out/'receipt.json').read_text())
    assert r['directed_dispositions']==4 and r['counts']=={'full:aligned':2,'plddt70:input_unavailable':2}
    with (out/'checkpoint_manifest.tsv').open() as f:
        checkpoints=list(csv.DictReader(f,delimiter='\t'))
    assert len(checkpoints)==4
    for checkpoint in checkpoints:assert sha(out/checkpoint['path'])==checkpoint['sha256']
    auditplan=dict(mode='reference',source_plan=str(pp),producer=dict(pid=99999999,created=0,cmdline=[]),output=str(root/'readback'),pins={str(pp):sha(pp)})
    ap=root/'audit-plan.json';write(ap,auditplan)
    subprocess.run([sys.executable,'scripts/readback_duplication_alignments.py','--plan',str(ap)],check=True)
    audit=json.loads((root/'readback/receipt.json').read_text())
    assert audit['directed_dispositions']==4 and audit['numerically_checked_alignments']==2
    assert audit['counts']==r['counts'] and audit['maximum_rmsd_rounding_error']<.00501
    # Preserve hashes but omit one disposition: exact full-grid verification must fail.
    checkpoint_manifest=out/'checkpoint_manifest.tsv'
    lines=checkpoint_manifest.read_text().splitlines();checkpoint_manifest.write_text('\n'.join(lines[:-1])+'\n')
    r['artifacts']['checkpoint_manifest.tsv']=sha(checkpoint_manifest);write(out/'receipt.json',r)
    auditplan['output']=str(root/'bad-readback');write(ap,auditplan)
    failure=subprocess.run([sys.executable,'scripts/readback_duplication_alignments.py','--plan',str(ap)],capture_output=True,text=True)
    assert failure.returncode!=0 and 'Incomplete disposition readback' in failure.stderr
print('Passed full numeric handoff with both native orders and excluded masks; missing disposition rejected despite updated hashes.')
