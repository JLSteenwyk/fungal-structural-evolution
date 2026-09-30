#!/usr/bin/env python3
"""Diagnose every reference alignment; preserve and quarantine RMSD discrepancies."""
import argparse
import csv
import hashlib
import json
import time
from collections import Counter
from functools import lru_cache
from pathlib import Path
import psutil
from duplication_alignment_numeric_readback import load_pdb
from duplication_alignment_numeric_diagnostic import check_alignment
from run_ortholog_pair_guide_comparison import sha
from expanded_reference_alignment_handoff import load_handoff, write_work_partition


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed readback plan')
        for path,h in plan['pins'].items():
            if sha(path)!=h:raise ValueError('Changed source: '+path)
    verify();dep=plan['producer']
    while psutil.pid_exists(dep['pid']):
        try:
            p=psutil.Process(dep['pid'])
            if abs(p.create_time()-dep['created'])>.01 or p.status()==psutil.STATUS_ZOMBIE:break
            if p.cmdline()!=dep['cmdline']:raise ValueError('Changed producer identity')
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();sp=Path(plan['source_plan']);sourceplan=json.loads(sp.read_text());root=Path(sourceplan['output'])
    rp=root/'receipt.json';receipt=json.loads(rp.read_text());receipt_hash=sha(rp)
    reference=plan['mode']=='reference'
    if plan['mode'] not in ['primary','reference']:raise ValueError('Unknown mode')
    status='complete_reference_alignment_dispositions_pending_readback' if reference else 'complete_duplication_alignment_dispositions_pending_readback'
    if receipt['status']!=status or receipt['plan_sha256']!=sha(sp):raise ValueError('Incomplete alignment source')
    for path,h in sourceplan['pins'].items():
        if sha(path)!=h:raise ValueError('Changed producer pin')
    if not reference:raise ValueError('Expanded reference mode required')
    verified_inputs,pairs,all_bindings,bundle,partition=load_handoff(sourceplan)
    if receipt['upstream_bindings']!=all_bindings:raise ValueError('Full upstream bindings differ')
    if (receipt['full_reference_pairs']!=len(partition)
            or receipt['existing_catalog_pairs_pending_reuse']!=len(partition)-len(pairs)):
        raise ValueError('Full reference work counts differ')
    work=root/'full_reference_work_partition.tsv'
    if sha(work)!=receipt['artifacts']['full_reference_work_partition.tsv']:raise ValueError('Changed full work partition')
    expected={}
    for row in pairs:
        ends=[(row['model_a'],int(row['version_a'])),(row['model_b'],int(row['version_b']))]
        pair=hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest()
        if pair!=row['pair_key'] or ends[0]==ends[1]:raise ValueError('Invalid source pair')
        for mask in ['full','plddt70']:
            for order,endpoints in enumerate([ends,ends[::-1]]):
                key=pair,mask,order
                if key in expected:raise ValueError('Duplicate source pair')
                expected[key]=endpoints
    if len(expected)!=receipt['directed_dispositions'] or len(pairs)!=receipt['distinct_model_pairs']:
        raise ValueError('Incomplete directed universe')
    specs=list(sourceplan['input_sources'].values()) if reference else [dict(inputs=sourceplan['inputs'],input_plan=sourceplan['input_plan'])]
    inputs={};bindings={}
    for spec in specs:
        folder=Path(spec['inputs']);ip=folder/'receipt.json';ir=json.loads(ip.read_text());manifest=folder/'inputs.jsonl'
        if ir['plan_sha256']!=sha(spec['input_plan']) or sha(manifest)!=ir['artifacts']['inputs.jsonl']:raise ValueError('Changed materialization')
        bindings[str(ip)]=sha(ip);bindings[str(manifest)]=sha(manifest)
        with manifest.open() as handle:
            for line in handle:
                row=json.loads(line);key=row['model_id'],row['version'],row['mask']
                if key in inputs:raise ValueError('Repeated input model/mask')
                inputs[key]=row
    if reference:
        if bindings!=receipt['input_bindings']:raise ValueError('Input bundle differs')
        mh=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest()
        if mh!=receipt['input_bundle_sha256']:raise ValueError('Input bundle digest differs')
    else:
        mh=sha(Path(sourceplan['inputs'])/'inputs.jsonl')
        if mh!=receipt['input_manifest_sha256'] or sha(Path(sourceplan['inputs'])/'receipt.json')!=receipt['input_receipt_sha256']:raise ValueError('Primary input binding differs')
    @lru_cache(maxsize=256)
    def coordinates(key):return load_pdb(inputs[key])
    cm=root/'checkpoint_manifest.tsv'
    if sha(cm)!=receipt['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed checkpoint manifest')
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    write_work_partition(out/'full_reference_work_partition.tsv',partition)
    if (out/'full_reference_work_partition.tsv').read_bytes()!=work.read_bytes():raise ValueError('Full work partition differs')
    seen=set();totals=Counter();max_error=0.;checked=0;discrepancies=[];rmsd_counts=Counter()
    table=out/'numeric_readback.tsv'
    fields=['pair_key','mask','order','rmsd_status','aligned_length','rmsd_recomputed','rmsd_native','rmsd_rounding_error','sequence_identity_exact','tm_left_native','tm_right_native','coverage_left','coverage_right','joint_plddt70_pairs','joint_plddt70_fraction']
    with cm.open() as handle,table.open('w') as resultfile:
        writer=csv.DictWriter(resultfile,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for item in csv.DictReader(handle,delimiter='\t'):
            path=root/item['path']
            if sha(path)!=item['sha256']:raise ValueError('Changed checkpoint')
            r=json.loads(path.read_text());key=r['pair_key'],r['mask'],r['order']
            if key in seen or key not in expected:raise ValueError('Unexpected/duplicate directed checkpoint')
            pair,mask,order=key
            if item['path']!=f'pairs/{pair[:2]}/{pair}-{mask}-{order}.json':raise ValueError('Checkpoint path differs')
            ends=expected[key];rows=[inputs[(*endpoint,mask)] for endpoint in ends]
            identities=[dict(model_id=e[0],version=e[1],status=row['status'],path=row.get('path'),sha256=row.get('sha256')) for e,row in zip(ends,rows)]
            ready=all(row['status']=='ready' for row in rows)
            command=[sourceplan['usalign'],*[row['path'] for row in rows],*sourceplan['options']] if ready else None
            if r['inputs']!=identities or r['command']!=command or r['plan_sha256']!=sha(sp) or r['input_manifest_sha256']!=mh or r['status']!=item['status']:raise ValueError('Checkpoint provenance differs')
            status=r['status']
            if status=='aligned':
                if not ready or r['returncode']!=0:raise ValueError('Invalid successful disposition')
                numeric=check_alignment(r,*[coordinates((*e,mask)) for e in ends])
                max_error=max(max_error,numeric['rmsd_rounding_error']);checked+=1
                rmsd_counts[numeric['rmsd_status']]+=1
                if numeric['rmsd_status']=='outside_printed_rounding':
                    discrepancies.append(dict(pair_key=pair,mask=mask,order=order,**numeric))
                writer.writerow(dict(pair_key=pair,mask=mask,order=order,**numeric))
            elif status=='input_unavailable':
                if ready or r['input_statuses']!=[row['status'] for row in rows]:raise ValueError('Invalid unavailable disposition')
            elif status=='native_error':
                if not ready or r['returncode']==0:raise ValueError('Invalid native error disposition')
            elif status=='parse_error':
                if not ready or r['returncode']!=0 or not r.get('error'):raise ValueError('Invalid parse error disposition')
            elif status=='timeout':
                if not ready or r['timeout_seconds']!=sourceplan['per_pair_timeout_seconds']:raise ValueError('Invalid timeout disposition')
            else:raise ValueError('Unknown disposition')
            seen.add(key);totals[mask+':'+status]+=1
            if len(seen)%10000==0:print('Checked dispositions',len(seen),'/',len(expected),flush=True)
    if seen!=set(expected) or dict(totals)!=receipt['counts']:raise ValueError('Incomplete disposition readback')
    verify()
    for path,h in all_bindings.items():
        if sha(path)!=h:raise ValueError('Changed upstream source during audit')
    for path,h in bindings.items():
        if sha(path)!=h:raise ValueError('Changed input during readback')
    if sha(rp)!=receipt_hash or sha(cm)!=receipt['artifacts']['checkpoint_manifest.tsv']:raise ValueError('Changed producer during readback')
    result=dict(status='complete_reference_alignment_rmsd_diagnostic_not_scientific_acceptance',scientific_eligibility=False,rmsd_status_counts=dict(rmsd_counts),rmsd_discrepancies=discrepancies,mode=plan['mode'],plan_sha256=ph,
                producer_receipt_sha256=receipt_hash,directed_dispositions=len(seen),numerically_checked_alignments=checked,
                counts=dict(totals),maximum_rmsd_rounding_error=max_error,full_reference_pairs=len(partition),existing_catalog_pairs_pending_reuse=len(partition)-len(pairs),artifacts={name:sha(out/name) for name in ['numeric_readback.tsv','full_reference_work_partition.tsv']},
                scope='Exact pair/order/mask grid, checkpoint hashes and provenance checked; all successful alignment mappings, '
                'sequence identities and least-squares RMSDs independently reconstructed from hashed PDBs. TM-scores checked '
                'against native text only, not reoptimized. Failed/excluded dispositions checked for internal consistency, '
                'not independently rerun or adjudicated. RMSD discrepancies are explicitly quarantined; strict original failure remains. Confidence uses rounded PDB values. No biological asymmetry inference.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
