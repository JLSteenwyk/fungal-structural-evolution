"""Bind the complete interval-pair grid to independently verified domain inputs."""
import csv,hashlib,json
from collections import Counter
from pathlib import Path
from run_ortholog_pair_guide_comparison import sha


def load_handoff(plan):
    inventory=Path(plan['inventory']);qr=json.loads((inventory/'receipt.json').read_text());folder=Path(plan['inputs']);rp=folder/'receipt.json';r=json.loads(rp.read_text());auditpath=Path(plan['input_readback']);audit=json.loads(auditpath.read_text())
    if audit['status']!='passed_full_duplication_domain_input_readback' or audit['plan_sha256']!=sha(plan['readback_plan']) or audit['producer_receipt_sha256']!=sha(rp):raise ValueError('Full domain input readback required')
    if r['status']!='complete_duplication_domain_input_materialization_pending_readback' or r['plan_sha256']!=sha(plan['input_plan']) or r['inventory_receipt_sha256']!=sha(inventory/'receipt.json'):raise ValueError('Domain input provenance differs')
    for name in ['intervals.jsonl','domain_pairs.tsv']:
        if sha(inventory/name)!=qr['artifacts'][name]:raise ValueError('Changed interval inventory')
    with (inventory/'intervals.jsonl').open() as f:intervals={r['interval_id']:r for line in f for r in [json.loads(line)]}
    manifest=folder/'inputs.jsonl';mh=sha(manifest)
    if mh!=r['artifacts']['inputs.jsonl']:raise ValueError('Changed domain manifest')
    inputs={};counts=Counter()
    with manifest.open() as f:
        for line in f:
            row=json.loads(line);key=row['interval_id'],row['mask']
            if key in inputs or key[0] not in intervals or key[1] not in ['full','plddt70']:raise ValueError('Unexpected interval input')
            expected=intervals[key[0]]
            if any(row[k]!=expected[k] for k in ['model_id','version','start','end','original_length','source_sha256']):raise ValueError('Domain input identity differs')
            if row['status'] not in ['ready','too_few_retained_residues','source_rejected']:raise ValueError('Unknown input disposition')
            inputs[key]={k:row[k] for k in ['model_id','version','start','end','status','path','sha256','sequence'] if k in row};counts[row['mask']+':'+row['status']]+=1
    if set(inputs)!={(iid,mask) for iid in intervals for mask in ['full','plddt70']} or dict(counts)!=r['counts'] or r['counts']!=audit['counts'] or len(inputs)!=r['input_dispositions'] or len(intervals)!=qr['unique_intervals']:raise ValueError('Domain input universe differs')
    with (inventory/'domain_pairs.tsv').open() as f:pairs=list(csv.DictReader(f,delimiter='\t'))
    seen=set()
    for row in pairs:
        ends=[row['interval_a'],row['interval_b']];key=hashlib.sha256(json.dumps(sorted(ends),separators=(',',':')).encode()).hexdigest()
        if key!=row['domain_pair_key'] or key in seen or ends[0]==ends[1] or not set(ends)<=intervals.keys():raise ValueError('Invalid domain pair')
        seen.add(key)
    if len(pairs)!=qr['unique_interval_pairs']:raise ValueError('Domain pair universe differs')
    bindings={str(path):sha(path) for path in [rp,manifest,auditpath]}
    bundle=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest()
    return inputs,pairs,bindings,bundle
