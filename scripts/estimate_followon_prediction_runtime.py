#!/usr/bin/env python3
"""Estimate short-queue runtime from frozen, length-stratified prediction receipts."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha, read_table


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','predictions','output']: parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--disposition-file',default='queue_disposition.tsv')
    parser.add_argument('--status-field',default='queue_disposition')
    parser.add_argument('--eligible-status',default='short_prediction_eligible')
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Use a new immutable runtime snapshot')
    inputs=checked_receipt(args.inputs)
    config_path=args.predictions/'config.json';config=json.loads(config_path.read_text());config_sha=sha(config_path)
    measurements=[]
    for path in sorted(args.predictions.glob('S*.json')):
        if path.name.endswith('.oom.json'): continue
        row=json.loads(path.read_text())
        if row['status']!='verified_prediction' or row['config_sha256']!=config_sha: raise ValueError('Mixed or unverified prediction configuration')
        seconds=float(row['inference_seconds'])
        if not math.isfinite(seconds) or seconds<=0: raise ValueError('Invalid runtime')
        measurements.append({'sequence_id':row['sequence_id'],'length':row['length'],'inference_seconds':seconds,
            'receipt_path':str(path),'receipt_sha256':sha(path)})
    pending=[r for r in read_table(args.inputs/args.disposition_file) if r[args.status_field]==args.eligible_status]
    summaries=[]
    for low,high in [(1,128),(129,256),(257,384),(385,512)]:
        times=[r['inference_seconds'] for r in measurements if low<=r['length']<=high]
        count=sum(low<=int(r['length'])<=high for r in pending)
        if count and len(times)<10: raise ValueError('Insufficient observed length-bin timing for projection')
        q=np.quantile(times,[.1,.5,.9]) if times else np.zeros(3)
        summaries.append({'length_min':low,'length_max':high,'observed_predictions':len(times),'pending_predictions':count,
            'observed_p10_seconds':float(q[0]),'observed_median_seconds':float(q[1]),'observed_p90_seconds':float(q[2]),
            'median_projected_gpu_hours':count*float(q[1])/3600,'p10_projected_gpu_hours':count*float(q[0])/3600,'p90_projected_gpu_hours':count*float(q[2])/3600})
    if sum(r['pending_predictions'] for r in summaries)!=len(pending): raise ValueError('Queue length scope differs')
    args.output.mkdir(parents=True)
    for name,rows in [('observed_prediction_timings.tsv',measurements),('length_bin_runtime.tsv',summaries)]:
        with (args.output/name).open('w') as handle:
            writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result={'status':'complete_short_queue_resource_projection','pending_short_predictions':len(pending),'observed_predictions':len(measurements),
        'median_projected_inference_hours':sum(r['median_projected_gpu_hours'] for r in summaries),
        'observed_quantile_projection_hours':[sum(r['p10_projected_gpu_hours'] for r in summaries),sum(r['p90_projected_gpu_hours'] for r in summaries)],
        'planning_wall_hours_with_50percent_overhead':1.5*sum(r['p90_projected_gpu_hours'] for r in summaries),
        'planned_gpu':'Existing authorized GPU after current work finishes; recheck availability before launch',
        'planning_vram_gb':24,'planning_output_gb':20,'cost':'Existing host; no paid resources',
        'disposition_selection':{'file':args.disposition_file,'field':args.status_field,'eligible_status':args.eligible_status},'input_receipt_sha256':sha(args.inputs/'receipt.json'),'prediction_config_sha256':config_sha,'script_sha256':sha(Path(__file__)),
        'interpretation':'Within-length-bin observed inference-time quantiles are planning scenarios, not uncertainty intervals. Output serialization, model loading and contention add overhead; 50 percent above the p90-based sum is an explicit allowance. Completed receipt timings are observational and not an independent artifact-quality audit. Applies only to current canonical <=512-residue configuration; no extrapolation to deferred long proteins or the whole atlas. This does not schedule or launch the next GPU run.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__': main()
