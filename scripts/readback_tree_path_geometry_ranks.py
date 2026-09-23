#!/usr/bin/env python3
"""Recompute every descriptive path/geometry rank row and equal-marker summary."""
import argparse
import csv
from collections import defaultdict
import hashlib
from itertools import groupby,product
import json
import math
from pathlib import Path
import shutil
from statistics import median
import time
import numpy as np
import psutil
from scipy.stats import spearmanr


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def rows(path):
    with Path(path).open() as f:yield from csv.DictReader(f,delimiter='\t')


def unique_table(path,keys):
    result={}
    for r in rows(path):
        key=tuple(r[k] for k in keys)
        if key in result:raise ValueError('Duplicate summary row')
        result[key]=r
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    a=ap.parse_args();plan=json.loads(a.plan.read_text());pins={str(a.plan):sha(a.plan),**plan['pins']}
    def verify():
        for p,h in pins.items():
            if sha(p)!=h:raise ValueError('Changed dependency: '+p)
    verify();out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False)
    if 'predecessor' in plan:
        (out/'state.json').write_text(json.dumps({'status':'waiting_for_full_benchmark'})+'\n');dep=plan['predecessor']
        while True:
            try:
                p=psutil.Process(dep['pid']);live=p.create_time()==dep['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
            except psutil.NoSuchProcess:live=False
            if not live:break
            time.sleep(20)
        cp=json.loads(Path(dep['plan']).read_text());rp=Path(cp['output'])/'receipt.json';r=json.loads(rp.read_text())
        if r['status']!='complete_all_cohort_descriptive_tree_path_geometry_benchmark' or r['plan_sha256']!=sha(dep['plan']):raise ValueError('Wrong or incomplete predecessor')
        for stage in r['stages']:
            if sha(stage['receipt'])!=stage['sha256']:raise ValueError('Changed benchmark stage')
        if plan['benchmark']!=cp['stages'][0]['output'] or plan['summary']!=cp['stages'][1]['output']:raise ValueError('Wrong source stage paths')
        pins[str(rp)]=sha(rp)
    verify()
    if psutil.virtual_memory().available<plan['resources']['minimum_available_memory_gib']*2**30 or shutil.disk_usage(out).free<plan['resources']['minimum_free_disk_gib']*2**30:raise ValueError('Insufficient resources')
    roots={k:Path(plan[k]) for k in ['benchmark','summary']};receipts={}
    for k,root in roots.items():
        rp=root/'receipt.json';pins[str(rp)]=sha(rp);receipts[k]=json.loads(rp.read_text())
        for name,h in receipts[k]['artifacts'].items():pins[str(root/name)]=h
    br=receipts['benchmark'];sr=receipts['summary']
    if sr['status']!='complete_descriptive_within_marker_rank_summary' or sr['source_receipt_sha256']!=sha(roots['benchmark']/'receipt.json') or br['status']!='complete_paired_site_tree_path_point_benchmark':raise ValueError('Invalid source binding')
    verify();(out/'state.json').write_text(json.dumps({'status':'recomputing_all_rank_rows'})+'\n')
    reported=unique_table(roots['summary']/'marker_rank_associations.tsv',['marker','path_metric','geometry_metric'])
    summaries=unique_table(roots['summary']/'equal_marker_summary.tsv',['path_metric','geometry_metric'])
    paths=['aa_tree_path_point','3di_af_tree_path_point','3di_af_empirical_tree_path_point','3di_llm_tree_path_point'];geoms=['ca_superposition_rmsd_angstrom','local_distance_mean_absolute_change_angstrom','pae10_local_mean_absolute_change_angstrom']
    seen=set();seen_rows=set();count=estimable=0;max_error=0.;aggregate=defaultdict(list)
    for marker,group in groupby(rows(roots['benchmark']/'path_geometry_points.tsv'),key=lambda r:r['marker']):
        if marker in seen:raise ValueError('Noncontiguous repeated marker group')
        seen.add(marker);data=list(group);count+=len(data)
        if len({(r['taxon_a'],r['taxon_b']) for r in data})!=len(data):raise ValueError('Duplicate source pair')
        flags={r['marker_review_status'] for r in data}
        if len(flags)!=1:raise ValueError('Inconsistent marker review labels')
        for path,geom in product(paths,geoms):
            key=marker,path,geom;row=reported[key];seen_rows.add(key);present=[r for r in data if r[geom]!='']
            x=[float(r[path]) for r in present];y=[float(r[geom]) for r in present]
            if any(not math.isfinite(v) for v in x+y):raise ValueError('Nonfinite source value')
            ok=len(x)>=3 and len(set(x))>1 and len(set(y))>1
            status='descriptive_estimable' if ok else 'insufficient_pairs_or_constant_values'
            if row['status']!=status or row['marker_review_status']!=next(iter(flags)) or int(row['pairs'])!=len(x) or int(row['omitted_geometry_pairs'])!=len(data)-len(x):raise ValueError('Wrong rank row scope/status')
            if ok:
                value=float(spearmanr(x,y).statistic);observed=float(row['rank_correlation']);error=abs(value-observed)
                if not math.isfinite(observed) or error>1e-12:raise ValueError('Rank correlation differs')
                aggregate[path,geom].append(value);estimable+=1;max_error=max(max_error,error)
            elif row['rank_correlation']!='':raise ValueError('Unestimable rank not missing')
    if seen_rows!=set(reported) or count!=br['accepted_pairs'] or len(seen)!=sr['markers'] or len(reported)!=sr['summary_rows']:raise ValueError('Incomplete rank grid')
    if set(summaries)!=set(product(paths,geoms)):raise ValueError('Incomplete equal-marker summary')
    for key,row in summaries.items():
        values=aggregate[key]
        if int(row['estimable_markers'])!=len(values) or int(row['negative_marker_correlations'])!=sum(v<0 for v in values):raise ValueError('Summary counts differ')
        if values:
            observed=float(row['median_marker_rank_correlation'])
            if not math.isfinite(observed) or abs(observed-median(values))>1e-12:raise ValueError('Summary median differs')
        elif row['median_marker_rank_correlation']!='':raise ValueError('Missing summary incorrectly populated')
    verify();result=dict(status='passed_all_descriptive_geometry_rank_rows_and_summaries',plan_sha256=sha(a.plan),markers=len(seen),accepted_pairs=count,rank_rows=len(reported),estimable_correlations=estimable,equal_marker_summary_rows=len(summaries),maximum_absolute_correlation_difference=max_error,source_hashes=pins,scope='Every rank-grid row, missing/constant case, denominator, review label and equal-marker median/count recomputed. SciPy Spearman used instead of producer rankdata/corrcoef call; shared ranking/numerical libraries, no p-values retained or independence/biological inference claimed.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');(out/'state.json').write_text(json.dumps({'status':result['status']})+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
