#!/usr/bin/env python3
"""Map a full audited local prediction queue after its conversion handoff completes."""
import argparse
from collections import Counter
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from advance_prediction_snapshot import ROOT, identity, sha


def checked(folder):
    row=json.loads((folder/'receipt.json').read_text())
    for name,digest in row['artifacts'].items():
        if sha(folder/name)!=digest:raise ValueError('Changed artifact: '+name)
    return row


def rows(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True)
    a=p.parse_args();c=json.loads(a.config.read_text());config_sha=sha(a.config)
    control=ROOT/c['control_output'];control.mkdir(parents=True,exist_ok=True)
    lock=(control/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (control/'launch.json').exists() or (control/'receipt.json').exists():
        raise FileExistsError('Inspect existing execution before recovery')
    def verify():
        if sha(a.config)!=config_sha:raise ValueError('Controller config changed')
        for name,digest in c['pinned_files'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Pinned dependency changed: '+name)
    verify();print('Waiting for completed conversion controller',c['predecessor_pid'],flush=True)
    while identity(c['predecessor_pid'])==c['predecessor_start_ticks']:time.sleep(20)
    verify()
    predecessor=ROOT/c['predecessor_control']/'receipt.json'
    handoff=json.loads(predecessor.read_text())
    if handoff['status']!='complete_prediction_snapshot_handoff' or handoff['config_sha256']!=c['predecessor_config_sha256']:
        raise ValueError('Prediction conversion did not complete cleanly')
    audit_path=ROOT/c['audit'];conversion_path=ROOT/c['conversion']
    audit=checked(audit_path);conversion=checked(conversion_path)
    if (sha(audit_path/'receipt.json')!=handoff['audit_receipt_sha256']
            or sha(conversion_path/'receipt.json')!=handoff['conversion_receipt_sha256']
            or conversion['source_audit_receipt_sha256']!=sha(audit_path/'receipt.json')
            or audit['partial_prediction_snapshot'] or audit['remaining_eligible']!=0
            or audit['predictions']!=c['expected_models'] or conversion['models']!=c['expected_models']):
        raise ValueError('Audit/conversion chain or model universe differs')
    inventory=conversion_path/'inventory.jsonl'
    records=[json.loads(line) for line in inventory.read_text().splitlines()]
    ids={r['sequence_id'] for r in rows(audit_path/'predictions.tsv')}
    if len(records)!=len(ids) or {r['record_id'] for r in records}!=ids:
        raise ValueError('Converted inventory identity grid differs')
    for record in records:
        if record['status']!='verified' or len(record['models'])!=1:raise ValueError('Unexpected converted model record')
        model=record['models'][0]
        if (model['sequence_sha256']!=record['record_id'][1:] or model['provider']!='local'
                or model['tool']!='ESMFold v1' or sha(ROOT/model['path'])!=model['sha256']):
            raise ValueError('Converted model identity/source/checksum differs')
    output=ROOT/c['mapping_output']
    if output.exists():raise FileExistsError('Use a fresh immutable mapping output')
    if shutil.disk_usage(ROOT).free<c['resource_plan']['output_headroom_bytes']:
        raise RuntimeError('Insufficient output headroom')
    command=[sys.executable,'scripts/map_marker_structures.py','--inventory',str(inventory),
             '--provider','local','--tool','ESMFold v1','--output',str(output)]
    launch={'command':command,'config_sha256':config_sha,'predecessor_receipt_sha256':sha(predecessor),
            'conversion_receipt_sha256':sha(conversion_path/'receipt.json'),'unix_time':time.time(),
            'resource_plan':c['resource_plan']}
    (control/'launch.json').write_text(json.dumps(launch,indent=2)+'\n')
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    with (control/'mapping.log').open('w') as log:
        subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    verify();mapping=checked(output)
    fields=['marker','taxon_id','protein_id','sequence_sha256']
    key=lambda r:tuple(r[k] for k in fields)
    expected=Counter(map(key,rows(audit_path/'taxon_links.tsv')))
    observed=Counter(map(key,rows(output/'marker_structure_links.tsv')))
    if expected!=observed or mapping['distinct_models']!=c['expected_models'] or mapping['source_inventory_sha256']!=sha(inventory):
        raise ValueError('Mapped model or originating marker-link universe differs')
    result={'status':'complete_marker_mapping_pending_residue_readback','config_sha256':config_sha,
            'mapping_receipt_sha256':sha(output/'receipt.json'),'models':mapping['distinct_models'],
            'marker_links':mapping['marker_proteins_linked'],'matrix_residue_links':mapping['matrix_residue_links'],
            'interpretation':'Full model/link identity handoff verified. Independent residue readback, PAE binding and native-feature qualification remain required before evolutionary use.'}
    (control/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
