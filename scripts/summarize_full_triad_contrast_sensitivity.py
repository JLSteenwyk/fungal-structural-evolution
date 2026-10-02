#!/usr/bin/env python3
"""Retain full contrast envelopes across masks, cores and every structural order."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil
from full_triad_contrast_sensitivity_sources import load_sources, selected_keys, MASKS, CORES, QUALIFIED
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def direction(n, expected, unique, low, high, tolerance):
    if n != expected: return 'unavailable'
    if not unique: return 'nonunique_fit'
    if low > tolerance: return 'positive'
    if high < -tolerance: return 'negative'
    if low >= -tolerance and high <= tolerance: return 'within_numerical_tolerance'
    return 'sign_uncertain'


def scenario(triad, leaves, mask, core, plan):
    keys = selected_keys(mask, core); chosen = [leaves[k] for k in keys]
    ranges = [g['metric_ranges']['rmsd_ar_minus_br'] for g in chosen]
    n = sum(g['available_orders'] for g in ranges); expected = 8 * len(keys)
    bounds = [(g['minimum'], g['maximum']) for g in ranges if g['available_orders']]
    low = min((a for a, b in bounds), default=None); high = max((b for a, b in bounds), default=None)
    unique = all(g['all_orders_unique'] for g in chosen)
    numerical = direction(n, expected, unique, low, high, plan['contrast_numerical_tolerance_angstrom'])
    strict = direction(n, expected, unique, low, high, 0.)
    screens = {}
    for s in plan['screens']:
        sid = s['id']; bits = [g['screens'][sid]['order_pass_bits'] for g in chosen]
        passed = all(b == 255 for b in bits)
        if passed: assert n == expected and unique
        screens[sid] = dict(leaf_order_pass_bits=bits, all_selected_orders_pass=passed,
                            qualified_direction=numerical if passed else 'excluded_by_quality')
    return dict(triad_id=triad['triad_id'], models=triad['models'], role_order=['a', 'b', 'reference'],
                mask_scenario=mask, core_scenario=core,
                source_order_group_keys=[[triad['triad_id'], m, c] for m, c in keys],
                expected_orders=expected, available_orders=n, all_selected_orders_unique=unique,
                minimum_contrast_angstrom=low, maximum_contrast_angstrom=high,
                contrast_span_angstrom=high-low if bounds else None,
                strict_zero_direction=strict, numerical_tolerance_direction=numerical, screens=screens)


def execute(plan, path, groups, triads, bindings, out):
    configuration = out/'configuration.json'; config = dict(plan_sha256=sha(path), source_hashes=bindings)
    if configuration.exists(): assert json.loads(configuration.read_text()) == config
    else: configuration.write_text(json.dumps(config, indent=2)+'\n')
    chunks = out/'checkpoints'; chunks.mkdir(exist_ok=True); expected_chunks=set(); artifacts={}
    counts=Counter(); qualified=Counter(); records=0; maximum=0.; reused=0; chunk=[]; chunk_number=0
    def save_chunk():
        nonlocal chunk_number, reused
        if not chunk: return
        target=chunks/f'{chunk_number:06d}.jsonl.gz'; expected_chunks.add(target.name)
        data=gzip.compress(('\n'.join(chunk)+'\n').encode(), compresslevel=1, mtime=0)
        if target.exists(): assert target.read_bytes()==data; reused+=1
        else:
            tmp=target.with_suffix('.tmp');tmp.write_bytes(data);tmp.replace(target)
        artifacts[str(target.relative_to(out))]=sha(target);chunk.clear();chunk_number+=1
    with gzip.open(groups,'rt') as source, gzip.open(out/'contrast_sensitivity.jsonl.gz','wt',compresslevel=1) as dest:
        for i,triad in enumerate(triads,1):
            leaves={}
            for m,c in itertools.product(MASKS[:2],CORES[:2]):
                line=next(source,None);assert line is not None;g=json.loads(line)
                assert [g['triad_id'],g['mask'],g['mapping_definition']]==[triad['triad_id'],m,c]
                assert g['models']==triad['models'] and g['role_order']==['a','b','reference']
                assert g['order_bits_layout']==[list(o) for o in itertools.product([0,1],repeat=3)]
                assert type(g['all_orders_unique']) is bool
                r=g['metric_ranges']['rmsd_ar_minus_br'];assert type(r['available_orders']) is int and 0<=r['available_orders']<=8
                assert (r['minimum'] is None)==(r['available_orders']==0)==(r['maximum'] is None)
                if r['available_orders']:assert math.isfinite(r['minimum']) and math.isfinite(r['maximum']) and r['minimum']<=r['maximum']
                for s in plan['screens']:
                    b=g['screens'][s['id']]['order_pass_bits'];assert type(b) is int and 0<=b<=255
                leaves[m,c]=g
            for m,c in itertools.product(MASKS,CORES):
                record=scenario(triad,leaves,m,c,plan);encoded=json.dumps(record,separators=(',',':'))
                dest.write(encoded+'\n');chunk.append(encoded);records+=1
                counts[f'{m}|{c}|'+record['numerical_tolerance_direction']]+=1
                if record['contrast_span_angstrom'] is not None:maximum=max(maximum,record['contrast_span_angstrom'])
                for sid,screen in record['screens'].items():qualified[f'{m}|{c}|{sid}|'+screen['qualified_direction']]+=1
            if i%plan['checkpoint_triads']==0:
                save_chunk();tmp=out/'state.tmp';tmp.write_text(json.dumps(dict(status='running_full_contrast_sensitivity',completed_triads=i,total_triads=len(triads)))+'\n');tmp.replace(out/'state.json')
                assert sum(p.stat().st_size for p in chunks.glob('*.jsonl.gz'))+dest.tell()<=plan['resources']['maximum_output_gib']*2**30
                assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            if i%1000==0:print('Full contrast sensitivity',i,'/',len(triads),flush=True)
        assert next(source,None) is None;save_chunk()
    assert {p.name for p in chunks.glob('*.jsonl.gz')}==expected_chunks
    fields=['mask_scenario','core_scenario','screen','qualified_direction','physical_triads','total_physical_triads']
    rows=0
    with (out/'contrast_sensitivity_counts.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for m,c,s in itertools.product(MASKS,CORES,plan['screens']):
            sid=s['id'];assert sum(qualified[f'{m}|{c}|{sid}|{d}'] for d in QUALIFIED)==len(triads)
            for d in QUALIFIED:
                writer.writerow(dict(zip(fields,[m,c,sid,d,qualified[f'{m}|{c}|{sid}|{d}'],len(triads)])));rows+=1
    summary=dict(measured_triads=len(triads),source_fit_rows=32*len(triads),source_order_groups=4*len(triads),
                 sensitivity_groups=records,screen_decisions=records*len(plan['screens']),summary_rows=rows,
                 direction_counts=dict(counts),qualified_direction_counts=dict(qualified),maximum_contrast_span=maximum)
    assert all(summary[k]==v for k,v in plan['expected'].items() if k in summary)
    bind(bindings,configuration);verify(bindings)
    for name in ['contrast_sensitivity.jsonl.gz','contrast_sensitivity_counts.tsv']:artifacts[name]=sha(out/name)
    result=dict(status='complete_full_triad_contrast_sensitivity_pending_independent_readback',**summary,
                plan_sha256=sha(path),source_hashes=bindings,artifacts=artifacts,checked_reused_chunks=reused,
                contrast_numerical_tolerance_angstrom=plan['contrast_numerical_tolerance_angstrom'],scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts','direction_counts','qualified_direction_counts']}),flush=True)
    return result


def run(path):
    plan=json.loads(Path(path).read_text());assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    raw,groups,triads,bindings=load_sources(plan,path);out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists()
        return execute(plan,path,groups,triads,bindings,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
