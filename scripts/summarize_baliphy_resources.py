#!/usr/bin/env python3
"""Summarize measured short-chain resources; never extrapolate convergence."""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import time
import psutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--wait-for-final',action='store_true')
    args=ap.parse_args();script_hash=sha(__file__)
    if args.wait_for_final:
        launch_path=Path('metadata/baliphy_sample_mapping_final_readback_launch_20260927.json')
        launch_hash=sha(launch_path);launch=json.loads(launch_path.read_text())
        while True:
            assert sha(__file__)==script_hash and sha(launch_path)==launch_hash
            state=dict(s.split('=',1) for s in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','MainPID','-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
            if state['ActiveState'] in ['inactive','failed']:
                assert state['ActiveState']=='inactive' and state['Result']=='success' and state['ExecMainStatus']=='0',state
                break
            assert state['ActiveState']=='active' and int(state['MainPID'])==launch['pid']
            try:
                p=psutil.Process(launch['pid'])
                assert p.create_time()==launch['created'] and p.cmdline()==launch['cmdline']
            except psutil.NoSuchProcess:
                time.sleep(1);continue
            time.sleep(30)
    audit_path=args.audit/'receipt.json';audit=json.loads(audit_path.read_text())
    for name,digest in audit['artifacts'].items():assert sha(args.audit/name)==digest
    eq=Path('results/ancestral/baliphy-input-equivalence-20260927-v2')
    eqr=json.loads((eq/'receipt.json').read_text())
    assert sha(eq/'configuration_mapping.json')==eqr['artifacts']['configuration_mapping.json']
    mapping={r['job_id']:r for r in json.loads((eq/'configuration_mapping.json').read_text())}
    dispositions=json.loads((args.audit/'dispositions.json').read_text())
    assert len(dispositions)==324 and {r['job_id'] for r in dispositions}==set(mapping)
    if args.wait_for_final:assert not any(r['status']=='pending_receipt' for r in dispositions)
    root=Path('results/ancestral/baliphy-sample-mapping-20260927-v1');rows=[];pins={str(audit_path):sha(audit_path),str(eq/'receipt.json'):sha(eq/'receipt.json')}
    for d in dispositions:
        row=dict(job_id=d['job_id'],group=mapping[d['job_id']]['group'],proteins=mapping[d['job_id']]['proteins'],audit_status=d['status'],producer_status='',elapsed_seconds='',peak_sampled_rss_bytes='',artifact_bytes='')
        if d['status']!='pending_receipt':
            path=root/d['job_id']/'receipt.json';digest=sha(path)
            assert digest==audit['pins'][str(path)];pins[str(path)]=digest
            r=json.loads(path.read_text())
            for name,h in r['artifacts'].items():assert sha(path.parent/name)==h
            row.update(producer_status=r['status'],elapsed_seconds=r.get('elapsed_seconds',''),peak_sampled_rss_bytes=r.get('peak_sampled_rss_bytes',''),artifact_bytes=sum((path.parent/name).stat().st_size for name in r['artifacts']))
        rows.append(row)
    groups=defaultdict(list)
    for row in rows:groups[row['group']].append(row)
    assert len(groups)==135
    group_rows=[]
    for group,members in sorted(groups.items()):
        measured=[r for r in members if r['elapsed_seconds']!='']
        group_rows.append(dict(group=group,labels=len(members),measured_attempts=len(measured),passed_samples=sum(r['audit_status']=='all_saved_samples_and_candidate_nodes_verified' for r in members),median_measured_seconds=statistics.median(r['elapsed_seconds'] for r in measured) if measured else '',maximum_sampled_rss_bytes=max((r['peak_sampled_rss_bytes'] for r in measured),default='')))
    args.output.mkdir(parents=True,exist_ok=False)
    for name,data in [('configuration_resources.tsv',rows),('effective_group_resources.tsv',group_rows)]:
        with (args.output/name).open('w') as h:
            w=csv.DictWriter(h,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    measured=[r for r in rows if r['elapsed_seconds']!='']
    result=dict(status='complete_resource_summary_of_frozen_audit',configurations=324,effective_input_groups=135,audit_dispositions=dict(Counter(r['audit_status'] for r in rows)),measured_attempts=len(measured),total_measured_worker_seconds=sum(r['elapsed_seconds'] for r in measured),maximum_sampled_rss_bytes=max((r['peak_sampled_rss_bytes'] for r in measured),default=None),pins=pins,script_sha256=script_hash,artifacts={p.name:sha(p) for p in args.output.iterdir()},scope='Elapsed time includes startup and only 20 iterations; sampled RSS can miss peaks. Capped attempts are censored, not successful runtime estimates. Alias groups retain all labels and do not represent independent chains. No production runtime, convergence, or posterior adequacy is inferred.')
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}))


if __name__=='__main__':main()
