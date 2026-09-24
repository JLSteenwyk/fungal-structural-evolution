#!/usr/bin/env python3
"""Summarize all candidate comparisons with conservative bidirectional metrics."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import defaultdict
from statistics import median


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['candidates','confidence','readback','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    sources=[Path(__file__)]
    for root in [args.candidates,args.confidence,args.readback]:sources.append(root/'receipt.json')
    candidate_receipt=json.loads((args.candidates/'receipt.json').read_text())
    confidence_receipt=json.loads((args.confidence/'receipt.json').read_text())
    readback_receipt=json.loads((args.readback/'receipt.json').read_text())
    if readback_receipt['status']!='passed_full_cross_clan_mapping_rmsd_identity_readback':raise ValueError('Full alignment readback required')
    for root,receipt in [(args.candidates,candidate_receipt),(args.confidence,confidence_receipt),(args.readback,readback_receipt)]:
        for name,h in receipt['artifacts'].items():
            if sha(root/name)!=h:raise ValueError('Artifact checksum differs')
            sources.append(root/name)
    pins={str(p):sha(p) for p in sources}
    candidates=rows(args.candidates/'candidates.tsv');members=rows(args.candidates/'candidate_members.tsv')
    intervals={r['interval_id']:r for r in rows(args.confidence/'interval_confidence.tsv')}
    directed={}
    for r in rows(args.readback/'alignment_readback.tsv'):
        key=(r['interval_left'],r['interval_right'])
        if key in directed:raise ValueError('Repeated alignment')
        directed[key]=r
    paired={}
    for a,b in sorted(directed):
        if a>=b:continue
        forward=directed[a,b];reverse=directed[b,a]
        minimum=lambda fields:min(float(r[f]) for r in [forward,reverse] for f in fields)
        paired[a,b]=dict(interval_a=a,interval_b=b,
            min_tm=minimum(['tm_left_native','tm_right_native']),
            min_coverage=minimum(['coverage_left','coverage_right']),
            min_joint_confidence=minimum(['joint_plddt70_fraction']),
            min_domain_confidence=min(float(intervals[i]['fraction_ca_plddt_ge70']) for i in [a,b]),
            max_rmsd=max(float(r['rmsd_recomputed']) for r in [forward,reverse]),
            max_tm_order_difference=max(abs(float(forward['tm_left_native'])-float(reverse['tm_right_native'])),abs(float(forward['tm_right_native'])-float(reverse['tm_left_native']))))
    if len(paired)*2!=len(directed):raise ValueError('Missing reverse comparison')
    keys=['boundary','representative','pfam_left','pfam_right'];groups=defaultdict(lambda:[set(),set()])
    for r in members:
        if r['model_exclusive_within_pair_cluster']=='1':groups[tuple(r[k] for k in keys)][r['side']=='right'].add(r['interval_id'])
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    pairtable=out/'bidirectional_pairs.tsv';summary=out/'candidate_confidence_summary.tsv'
    with pairtable.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(next(iter(paired.values()))),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(paired.values())
    fields=keys+['confidence_fraction','all_interval_pairs','retained_interval_pairs','screen_interval_pairs',
                 'median_min_tm','median_min_coverage','median_max_rmsd','max_tm_order_difference']
    seen=set();summaries=[]
    for candidate in candidates:
        key=tuple(candidate[k] for k in keys)
        if key in seen:raise ValueError('Repeated candidate')
        seen.add(key);left,right=groups[key]
        ids={tuple(sorted((a,b))) for a in left for b in right}
        values=[paired[i] for i in sorted(ids)]
        for threshold in [0,.5,.8,.9]:
            kept=[v for v in values if min(v['min_domain_confidence'],v['min_joint_confidence'])>=threshold]
            screened=[v for v in kept if v['min_tm']>=.5 and v['min_coverage']>=.8]
            row=dict(zip(keys,key));row.update(confidence_fraction=threshold,all_interval_pairs=len(values),
                retained_interval_pairs=len(kept),screen_interval_pairs=len(screened))
            for field,metric in [('median_min_tm','min_tm'),('median_min_coverage','min_coverage'),('median_max_rmsd','max_rmsd')]:
                row[field]=median(v[metric] for v in kept) if kept else ''
            row['max_tm_order_difference']=max(v['max_tm_order_difference'] for v in kept) if kept else ''
            summaries.append(row)
    with summary.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summaries)
    for p,h in pins.items():
        if sha(p)!=h:raise ValueError('Source changed')
    result=dict(status='complete_descriptive_cross_clan_confidence_summary',source_hashes=pins,
                unordered_pairs=len(paired),candidate_entries=len(candidates),summary_rows=len(summaries),
                screen_entries_by_threshold={str(t):sum(r['confidence_fraction']==t and r['screen_interval_pairs']>0 for r in summaries) for t in [0,.5,.8,.9]},
                artifacts={p.name:sha(p) for p in [pairtable,summary]},
                scope='Descriptive only. Each pair uses the minimum score/coverage/confidence across both input orders and length normalizations. '
                'Confidence fraction applies to both full domains and both matched-residue alignments. '
                'Screen requires minimum TM-score >=0.5 and coverage >=0.8; these are prioritization choices, not proof of homology. '
                'Boundary views, source models and pair observations overlap; no independent tests, PAE filtering, function or evolutionary claims.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
