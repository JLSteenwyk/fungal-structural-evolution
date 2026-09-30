#!/usr/bin/env python3
"""Native two-source handoff, primary-pair reuse and corruption fixture."""
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
from expanded_reference_alignment_handoff import load_handoff
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
    write(ledger,dict(status='passed_full_reference_comparison_ledger_and_native_model_readback',additional_models=1,
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
                         status='ready',path=str(pdb),sha256=sha(pdb),sequence=seq,original_positions=list(range(1,len(seq)+1)),retained_residues=len(seq),original_length=len(seq))
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
        proof=root/(label+'-proof.json')
        status='passed_full_duplication_alignment_input_readback' if label=='primary' else 'passed_full_reference_alignment_input_readback'
        field='source_receipt_sha256' if label=='primary' else 'producer_receipt_sha256'
        write(proof,dict(status=status,**{field:sha(folder/'receipt.json')}))
        sources[label]=dict(inputs=str(folder),input_plan=str(ip),readback=str(proof))
    reuse=root/'reuse';reuse.mkdir()
    planning=[]
    for a,b,disp in [('A','B','matching_catalog_sources_pending_input_and_result_checks'),('A','C','new_pair_requires_alignment')]:
        row=pair(a,b,'ignored');row.pop('work_disposition');row.update(new_sources=json.dumps(['reference_expanded']),matching_old_sources=json.dumps(['primary_expanded_completed'] if a=='A' and b=='B' else []),changed_old_sources='[]',disposition=disp);planning.append(row)
    table(reuse/'pair_reuse_candidates.tsv',planning)
    write(reuse/'receipt.json',dict(status='complete_full_pair_union_catalog_reuse_screen_not_authorization',artifacts={'pair_reuse_candidates.tsv':sha(reuse/'pair_reuse_candidates.tsv')},source_hashes={str(inventory/n):sha(inventory/n) for n in ['receipt.json','models.jsonl','model_pairs.tsv']},per_new_source_counts={'reference_expanded:matching_catalog_sources_pending_input_and_result_checks':1,'reference_expanded:new_pair_requires_alignment':1}))
    reuseproof=root/'reuse-proof.json';write(reuseproof,dict(status='passed_full_pair_union_reuse_candidate_sql_readback',producer_receipt_sha256=sha(reuse/'receipt.json')))
    plan=dict(reuse_candidates=str(reuse),reuse_readback=str(reuseproof),inventory=str(inventory),base_queue=str(base),inventory_readback=str(ledger),input_sources=sources,
              output=str(root/'out'),workers=2,minimum_free_disk_gib=0,
              usalign='/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/US-align/USalign',
              options=['-mol','prot','-mm','0','-outfmt','0','-ter','2'],per_pair_timeout_seconds=10,pins={})
    inputs,pairs,bindings,bundle,partition=load_handoff(plan)
    assert len(partition)==2 and partition[0]['measurement_disposition'] in ['native_measurement_required','pending_actual_input_result_and_numeric_reuse_checks']
    assert len(inputs)==6 and len(pairs)==1 and pairs[0]['model_b']=='C'
    pp=root/'plan.json';write(pp,plan)
    subprocess.run([sys.executable,'scripts/run_expanded_reference_alignments.py','--plan',str(pp)],check=True)
    out=Path(plan['output']);r=json.loads((out/'receipt.json').read_text())
    assert r['full_reference_pairs']==2 and r['existing_catalog_pairs_pending_reuse']==1
    assert r['directed_dispositions']==4 and r['counts']=={'full:aligned':2,'plddt70:input_unavailable':2}
    with (out/'checkpoint_manifest.tsv').open() as f:
        checkpoints=list(csv.DictReader(f,delimiter='\t'))
    assert len(checkpoints)==4
    for checkpoint in checkpoints:assert sha(out/checkpoint['path'])==checkpoint['sha256']
    # Independent native-text/residue/least-squares diagnostic covers every new job.
    producer=dict(pid=99999999,created=0,cmdline=[])
    dp=root/'diagnostic-plan.json'
    write(dp,dict(source_plan=str(pp),producer=producer,mode='reference',pins={},output=str(root/'diagnostic')))
    subprocess.run([sys.executable,'scripts/diagnose_expanded_reference_alignments.py','--plan',str(dp)],check=True)
    dr=json.loads((root/'diagnostic/receipt.json').read_text())
    assert dr['directed_dispositions']==4 and dr['numerically_checked_alignments']==2
    assert dr['full_reference_pairs']==2 and dr['existing_catalog_pairs_pending_reuse']==1
    # Rehash a false checkpoint provenance field: the full diagnostic must reject it.
    checkpoint=out/checkpoints[0]['path'];original_checkpoint=checkpoint.read_bytes()
    altered=json.loads(original_checkpoint);altered['order']=7;write(checkpoint,altered)
    checkpoints[0]['sha256']=sha(checkpoint);table(out/'checkpoint_manifest.tsv',checkpoints)
    r['artifacts']['checkpoint_manifest.tsv']=sha(out/'checkpoint_manifest.tsv');write(out/'receipt.json',r)
    bad=root/'diagnostic-bad-plan.json';write(bad,dict(source_plan=str(pp),producer=producer,mode='reference',pins={},output=str(root/'diagnostic-bad')))
    result=subprocess.run([sys.executable,'scripts/diagnose_expanded_reference_alignments.py','--plan',str(bad)],capture_output=True)
    assert result.returncode!=0 and not (root/'diagnostic-bad/receipt.json').exists()
    # Rehash an incomplete planning union; endpoint/source verification must reject.
    original_planning=(reuse/'pair_reuse_candidates.tsv').read_bytes();original_reuse=(reuse/'receipt.json').read_bytes()
    table(reuse/'pair_reuse_candidates.tsv',planning[1:])
    rr=json.loads(original_reuse);rr['artifacts']['pair_reuse_candidates.tsv']=sha(reuse/'pair_reuse_candidates.tsv');write(reuse/'receipt.json',rr)
    write(reuseproof,dict(status='passed_full_pair_union_reuse_candidate_sql_readback',producer_receipt_sha256=sha(reuse/'receipt.json')))
    try:load_handoff(plan)
    except ValueError:pass
    else:raise AssertionError('Incomplete full work partition accepted')
    (reuse/'pair_reuse_candidates.tsv').write_bytes(original_planning);(reuse/'receipt.json').write_bytes(original_reuse)
    write(reuseproof,dict(status='passed_full_pair_union_reuse_candidate_sql_readback',producer_receipt_sha256=sha(reuse/'receipt.json')))
    # Rehash a modified manifest to ensure source provenance is checked beyond its file hash.
    manifest=root/'additional/inputs.jsonl';rows=[json.loads(l) for l in manifest.read_text().splitlines()]
    rows[0]['source_sha256']='wrong';manifest.write_text(''.join(json.dumps(r)+'\n' for r in rows))
    rp=root/'additional/receipt.json';r=json.loads(rp.read_text());r['artifacts']['inputs.jsonl']=sha(manifest);write(rp,r)
    try:load_handoff(plan)
    except ValueError as e:assert 'source hash' in str(e)
    else:raise AssertionError('Changed source provenance accepted')
print('Passed full source partition, verified two-source handoff, both masks/orders, native checkpoint/coordinate diagnostic, and rehashed false checkpoint/incomplete partition/source rejection.')
