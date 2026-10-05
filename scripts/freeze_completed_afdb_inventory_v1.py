#!/usr/bin/env python3
"""Freeze and validate the completed full retrieval log without restarting retrieval."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import fcntl
import json
import math
from pathlib import Path
import re
import subprocess

from ancestral_chain_attempt import sha
from freeze_append_only_inventory import freeze
from reference_measurement_union_sources import bind,verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.receipt.exists();pins={}
    plan=Path('metadata/recovery_20260925_structure_retrieval_plan.json')
    launch=Path('metadata/recovery_20260925_structure_retrieval_launch.json')
    pp,lp=[json.loads(q.read_text()) for q in [plan,launch]]
    final=Path(pp['output'])/'receipt.json';outcome=json.loads(final.read_text())
    assert outcome['status']=='complete_snapshot_retrieval_queue' and not outcome['stopped_for_disk_headroom']
    assert outcome['plan_sha256']==sha(plan) and outcome['last_log_line']=='Snapshot queue completed 1589879'
    state=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',lp['unit'],
        '-p','LoadState','-p','ActiveState','-p','Result','-p','ExecMainPID','-p','ExecMainCode','-p','ExecMainStatus',
        '-p','ExecMainStartTimestamp','-p','ExecMainExitTimestamp'],text=True).splitlines())
    # This older service was not collected: require its retained actual ExecMain identity/status.
    assert state['LoadState']=='loaded' and state['ActiveState']=='inactive'
    assert state['Result']=='success' and int(state['ExecMainPID'])==lp['pid']==863281
    assert state['ExecMainCode']=='1' and state['ExecMainStatus']=='0'
    assert state['ExecMainStartTimestamp'] and state['ExecMainExitTimestamp']
    for q in [plan,launch,final,Path(__file__),Path('scripts/freeze_append_only_inventory.py'),
              Path('scripts/retrieve_matched_models.py'),Path('scripts/recover_structure_retrieval.py')]:bind(pins,q)
    log=Path('data/raw/afdb_models.jsonl');lock=Path('data/structures/afdb/.retrieval.lock')
    with lock.open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        frozen=freeze(log,a.output)
        assert frozen['incomplete_tail_bytes_excluded']==0
        assert frozen['frozen_bytes']==frozen['source_initial_bytes']==frozen['source_observed_final_bytes']
    statuses=Counter();records=0;models=0;target_models=0
    with (a.output/'inventory.jsonl').open() as f:
        for line in f:
            row=json.loads(line);assert re.fullmatch(r'[A-Za-z0-9]+',row['uniprot_accession'])
            assert row['status'] in ['verified','no_exact_full_length_model','error'];statuses[row['status']]+=1;records+=1
            if row['status']=='verified':
                assert row['models']
                for model in row['models']:
                    models+=1
                    if (model.get('provider'),model.get('tool'))!=('GDM','AlphaFold Monomer v2.0 pipeline'):continue
                    target_models+=1
                    assert type(model['length']) is int and model['length']>0
                    assert type(model['version']) is int and model['version']>0
                    for key in ['sequence_sha256','sha256']:assert re.fullmatch(r'[0-9a-f]{64}',model[key])
                    assert math.isfinite(model['mean_ca_plddt']) and 0<=model['mean_ca_plddt']<=100
                    assert math.isfinite(model['fraction_ca_plddt_below50']) and 0<=model['fraction_ca_plddt_below50']<=1
                    assert model['path'].startswith('data/structures/afdb/') and '..' not in Path(model['path']).parts
    assert records==frozen['complete_lines']
    for q in [a.output/'receipt.json',a.output/'inventory.jsonl']:bind(pins,q)
    verify(pins)
    result=dict(status='passed_completed_afdb_inventory_freeze_and_record_contracts',checked_utc=datetime.now(timezone.utc).isoformat(),
        inventory=str(a.output/'inventory.jsonl'),inventory_sha256=frozen['inventory_sha256'],inventory_bytes=frozen['frozen_bytes'],
        complete_records=records,append_record_status_counts=dict(statuses),append_model_occurrences=models,
        append_gdm_model_occurrences=target_models,completed_original_queued_accessions=1589879,
        retained_original_service_state=state,source_hashes=pins,scientific_eligibility=False,new_downloads=0,gpu=False,
        scope='Original non-collected retrieval service actual retained ExecMain PID/code/status plus full terminal queue '
              'receipt checked. Exclusive original retrieval lock acquired; entire completed log copied and source/target '
              'prefix rehashed with zero excluded tail. All record statuses and target-source ranking/identity/length '
              'contracts checked. Append-record/model occurrences are not distinct current accessions, coverage or '
              'independent proteins. Full latest-status/model ranking, representative sequence linkage, coordinate-byte '
              'hashing and independent full-catalog readback remain required. No retrieval restart or structure prediction.')
    with a.receipt.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
