#!/usr/bin/env python3
"""Inventory all four-chain diagnostic groups and uncompressed planning scenarios."""
import argparse
from collections import defaultdict,Counter
import json
from pathlib import Path
from Bio import SeqIO
from ancestral_chain_attempt import sha,write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--producer',type=Path,required=True);ap.add_argument('--extraction',type=Path,required=True);ap.add_argument('--benchmark',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    p=json.loads(args.producer.read_text());e=json.loads(args.extraction.read_text());b=json.loads(args.benchmark.read_text())
    for plan in [p,e,b]:
        for path,h in plan['pins'].items():assert sha(path)==h,path
    jobs=json.loads(Path(p['jobs']).read_text());groups=defaultdict(list)
    for j in jobs:groups[j['config']['model_input_identity']].append(j)
    assert len(jobs)==1620 and len(groups)==405
    rows=[];lengths={}
    for key,js in sorted(groups.items()):
        assert len(js)==len({j['chain']['seed'] for j in js})==4
        assert len({j['chain']['alignment_sha256'] for j in js})==1
        chain=js[0]['chain'];path=chain['alignment'];h=chain['alignment_sha256']
        if path not in lengths:
            assert sha(path)==h
            records=list(SeqIO.parse(path,'fasta'));assert len(records)==len({r.id for r in records})
            lengths[path]=sum(len(str(r.seq).replace('-','')) for r in records)
        available=[]
        for j in js:
            dp=Path(e['output'])/j['chain']['chain_id']/'disposition.json'
            available.append(dp.exists() and json.loads(dp.read_text())['status']=='state_trace_complete_not_posterior_qualification')
        rows.append(dict(group=key,coordinates=4*lengths[path],chain_ids=[j['chain']['chain_id'] for j in js],extracted_chains=sum(available),candidate_ready=all(available)))
    coordinates=sum(r['coordinates'] for r in rows)
    per_draw={draws:max(r['maximum_seconds'] for r in b['rows'] if r['draws_per_chain']==draws) for draws in [50,75]}
    # Every coordinate distinct is an exact upper limit on pattern count, but
    # benchmark timing and serialized-size extrapolations are not hard bounds.
    result=dict(status='full_categorical_group_resource_inventory',groups=rows,group_count=len(rows),chains=len(jobs),
        coordinates_per_cutoff=coordinates,maximum_coordinates_per_group=max(r['coordinates'] for r in rows),
        candidate_ready_groups=sum(r['candidate_ready'] for r in rows),extracted_chain_count=sum(r['extracted_chains'] for r in rows),
        no_compression_pattern_count_two_cutoffs=2*coordinates,
        synthetic_maximum_timing_cpu_hours_no_compression=coordinates*sum(per_draw.values())/3600,
        synthetic_maximum_raw_report_gib_no_compression=2*coordinates*max(r['raw_bytes_per_pattern'] for r in b['rows'])/1024**3,
        pins={str(path):sha(path) for path in [args.producer,args.extraction,args.benchmark,Path(p['jobs']),Path(__file__)]},
        scope='All 405 groups. Readiness is file/status screening only; full provenance checks remain mandatory before diagnostics. Every coordinate counted, without extrapolating single-chain compression to quartets. Synthetic per-pattern performance is a planning scenario, not a guaranteed upper runtime or disk bound.')
    assert not args.output.exists();write_json(args.output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['groups','pins']},indent=2))


if __name__=='__main__':main()
