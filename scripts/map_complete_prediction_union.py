#!/usr/bin/env python3
"""Map an independently verified model union and audit all retained residues."""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from audit_busco_gene_copies import ROOT, sha, read_table


def checked(folder):
    receipt=json.loads((folder/'receipt.json').read_text())
    for name,digest in receipt['artifacts'].items():
        if sha(folder/name)!=digest:
            raise ValueError('Changed artifact: '+str(folder/name))
    return receipt


def table(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args()
    plan=json.loads(a.plan.read_text());plan_sha=sha(a.plan)
    def verify():
        if sha(a.plan)!=plan_sha:
            raise ValueError('Plan changed')
        for name,digest in plan['pins'].items():
            if sha(ROOT/name)!=digest:
                raise ValueError('Changed dependency: '+name)
    verify()
    inventory=ROOT/plan['inventory'];audit=ROOT/plan['inventory_readback']
    ir=checked(inventory);ar=checked(audit)
    if (ir['status']!='complete_disjoint_prediction_inventory_union_with_batch_provenance'
            or ar['status']!='passed_complete_disjoint_prediction_inventory_readback'
            or ar['producer_receipt_sha256']!=sha(inventory/'receipt.json')
            or ar['models']!=plan['models'] or ir['models']!=plan['models']
            or ar['matched_marker_protein_links']!=plan['marker_links']):
        raise ValueError('Verified complete inventory required')
    control,combined,mapping,readback=[ROOT/plan[k] for k in ['control','union_audit','mapping','readback']]
    if any(p.exists() for p in [control,combined,mapping,readback]):
        raise FileExistsError('Use fresh immutable outputs')
    if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient disk')
    control.mkdir(parents=True)
    state=dict(status='binding_prediction_audits',plan_sha256=plan_sha,started_unix=time.time())
    def save():
        state['updated_unix']=time.time()
        p=control/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(control/'state.json')
    save()
    records={row['record_id']:row for row in (json.loads(line) for line in (inventory/'inventory.jsonl').open())}
    predictions=[];links=[];seen=set();sources=[]
    for source in plan['audits']:
        folder=ROOT/source['path'];receipt=checked(folder)
        if receipt['status']!=source['status'] or receipt['predictions']!=source['models']:
            raise ValueError('Prediction audit scope/status differs')
        if receipt['status']=='complete_artifact_readback' and (receipt['partial_prediction_snapshot'] or receipt['remaining_eligible']):
            raise ValueError('Incomplete individual prediction audit')
        rows=read_table(folder/'predictions.tsv')
        if len(rows)!=source['models']:
            raise ValueError('Audit row count differs')
        local=set()
        for row in rows:
            sid=row['sequence_id']
            if sid in seen or sid not in records:
                raise ValueError('Overlapping or unexpected prediction')
            model=records[sid]['models'][0]
            if row['prediction_receipt_sha256']!=model['prediction_receipt_sha256'] or int(row['length'])!=model['length']:
                raise ValueError('Prediction audit/model binding differs')
            seen.add(sid);local.add(sid)
        local_links=read_table(folder/'taxon_links.tsv')
        if any(row['sequence_id'] not in local or row['sequence_id']!='S'+row['sequence_sha256'] for row in local_links):
            raise ValueError('Unexpected source links')
        predictions.extend(rows);links.extend(local_links)
        sources.append(dict(path=source['path'],receipt_sha256=sha(folder/'receipt.json'),
                            original_receipt=receipt,link_columns=list(local_links[0])))
    if seen!=set(records):
        raise ValueError('Missing audited models')
    expected=read_table(audit/'expected_global_links.tsv')
    key=lambda r:tuple(r[k] for k in ['marker','taxon_id','protein_id','sequence_sha256'])
    if Counter(map(key,links))-Counter(map(key,expected)):
        raise ValueError('Originating links outside exact global inventory')
    combined.mkdir(parents=True)
    table(combined/'predictions.tsv',sorted(predictions,key=lambda r:r['sequence_id']))
    table(combined/'taxon_links.tsv',links)
    (combined/'source_audits.json').write_text(json.dumps(sources,indent=2)+'\n')
    cr=dict(status='complete_disjoint_prediction_audit_union',predictions=len(seen),taxon_marker_links=len(links),
            artifacts={p.name:sha(p) for p in combined.iterdir()},
            scope='Previously audited sources bound to the complete model inventory. Original audit statuses and schemas retained; absent source columns are blank, not imputed. No new prediction or shared configuration claimed.')
    (combined/'receipt.json').write_text(json.dumps(cr,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    commands=[('mapping',[sys.executable,'scripts/map_marker_structures.py','--inventory',str(inventory/'inventory.jsonl'),
                          '--provider','local','--tool','ESMFold v1','--output',str(mapping)]),
              ('residue_readback',[sys.executable,'scripts/audit_local_residue_mapping_global_links.py','--mapping',str(mapping),
                                   '--prediction-audit',str(combined),'--expected-links',str(audit/'expected_global_links.tsv'),'--output',str(readback)])]
    for stage,command in commands:
        verify()
        if shutil.disk_usage(ROOT).free<plan['resources']['minimum_free_disk_gib']*2**30:
            raise ValueError('Insufficient disk')
        state.update(status='running_stage',stage=stage,command=command);save()
        with (control/(stage+'.log')).open('w') as log:
            subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    mr=checked(mapping);rr=checked(readback)
    if (mr['distinct_models']!=plan['models'] or mr['source_inventory_sha256']!=sha(inventory/'inventory.jsonl')
            or rr['status']!='passed_full_local_residue_mapping_readback' or rr['models']!=plan['models']
            or rr['marker_links']!=plan['marker_links'] or rr['mapping_receipt_sha256']!=sha(mapping/'receipt.json')):
        raise ValueError('Final mapping/readback scope differs')
    verify()
    result=dict(status='complete_full_prediction_union_mapping_and_residue_readback',models=rr['models'],marker_links=rr['marker_links'],
                matrix_residue_links=rr['matrix_residue_links'],mapping_receipt_sha256=sha(mapping/'receipt.json'),
                readback_receipt_sha256=sha(readback/'receipt.json'),union_audit_receipt_sha256=sha(combined/'receipt.json'),plan_sha256=plan_sha,
                scope='Every completed model mapped with all exact global links and independent residue readback. Batch configurations remain explicit. Qualified encoding union and paired evolutionary inference remain pending.')
    (control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    state.update(status=result['status']);save()


if __name__=='__main__':
    main()
