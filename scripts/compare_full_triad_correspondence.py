#!/usr/bin/env python3
"""Retain every sequence/structural order pair and summarize dependent alternatives."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import json
import math
from pathlib import Path
import shutil

from full_triad_correspondence_comparison_sources import load_sources,blocks,NUMERIC_FIELDS,SEQUENCE_NUMERIC_FIELDS,METHODS,PERMUTATIONS,DEFINITIONS,pair_fields
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def ranges(rows,fields,availability):
    result={}
    for field in fields:
        numbers=[float(r[field]) for r in rows if r[field]!=''];assert all(math.isfinite(v) for v in numbers)
        result[field]=dict(**{availability:len(numbers)},minimum=min(numbers) if numbers else None,
            maximum=max(numbers) if numbers else None,span=max(numbers)-min(numbers) if numbers else None,
            mean=math.fsum(numbers)/len(numbers) if numbers else None)
    return result


def sequence_summary(triad,mask,method,rows,plan):
    statuses=[r['fit_status'] for r in rows];unique=all(s=='computed_unique_at_numeric_tolerance' for s in statuses)
    metrics=ranges(rows,SEQUENCE_NUMERIC_FIELDS,'available_orders');contrast=metrics['rmsd_ar_minus_br']
    direction='uncomputed_or_excluded_order' if contrast['available_orders']!=6 else 'nonunique_order_fit' if not unique else 'positive' if contrast['minimum']>0 else 'negative' if contrast['maximum']<0 else 'includes_zero'
    screens={}
    for s in plan['screens']:
        sid=s['id'];values=[[int(r[sid+suffix]) for r in rows] for suffix in ['_pass','_core_pass','_three_pair_pass']]
        assert all(v in [0,1] for group in values for v in group) and all(a==(b and c) for a,b,c in zip(*values))
        bits=[sum(value<<i for i,value in enumerate(flags)) for flags in values]
        screens[sid]=dict(order_pass_bits=bits[0],core_order_pass_bits=bits[1],inherited_order_pass_bits=bits[2],
            passing_orders=sum(values[0]),all_orders_pass=bits[0]==63,any_order_pass=bits[0]!=0,
            order_exclusions=[r[sid+'_exclusions'].split(';') if r[sid+'_exclusions'] else [] for r in rows])
        if bits[0]==63:assert unique
    return dict(triad_id=triad['triad_id'],models=triad['models'],mask=mask,sequence_method=method,order_layout=PERMUTATIONS,
        order_fit_statuses=statuses,order_native_statuses=[r['sequence_native_status'] for r in rows],
        order_native_payload_hashes=[r['native_payload_sha256'] for r in rows],order_triplet_hashes=[r['triples_sha256'] for r in rows],
        distinct_order_triplet_sets=len({r['triples_sha256'] for r in rows}),all_orders_unique=unique,
        strict_all_orders_contrast_direction=direction,metric_ranges=metrics,screens=screens)


def pair(triad,mask,method,permutation,definition,index,s,t,st,tt,plan):
    left=set(st);right=set(tt);intersection=len(left&right);union=len(left|right)
    unique=s['fit_status']==t['fit_status']=='computed_unique_at_numeric_tolerance'
    sign=lambda r:1 if float(r['rmsd_ar_minus_br'])>0 else -1 if float(r['rmsd_ar_minus_br'])<0 else 0
    row=dict(triad_id=triad['triad_id'],mask=mask,sequence_method=method,sequence_permutation=permutation,mapping_definition=definition,
        structural_order=f'{index:03b}',sequence_fit_status=s['fit_status'],structural_fit_status=t['fit_status'],
        sequence_triples_sha256=s['triples_sha256'],structural_triples_sha256=t['triples_sha256'],
        sequence_common_residues=len(st),structural_common_residues=len(tt),intersection_triples=intersection,union_triples=union,
        triple_jaccard=intersection/union if union else '',identical_correspondence=int(left==right),both_unique_numeric=int(unique),
        strict_numeric_contrast_sign_agreement=('same' if sign(s)==sign(t) else 'different') if unique else 'unavailable')
    for field in NUMERIC_FIELDS:row['sequence_minus_structural_'+field]=float(s[field])-float(t[field]) if s[field]!='' and t[field]!='' else ''
    for screen in plan['screens']:
        sid=screen['id'];a=int(s[sid+'_pass']);b=int(t[sid+'_pass']);assert a in [0,1] and b in [0,1]
        row.update({sid+'_sequence_pass':a,sid+'_structural_pass':b,sid+'_both_pass':a*b})
    return row


def joint_summary(triad,mask,method,definition,rows,plan):
    names=['sequence_minus_structural_'+f for f in NUMERIC_FIELDS]
    metrics=ranges(rows,names,'available_pairs');contrast=metrics['sequence_minus_structural_rmsd_ar_minus_br']
    unique=sum(r['both_unique_numeric'] for r in rows)
    direction='uncomputed_or_nonunique_pair' if unique!=48 else 'positive' if contrast['minimum']>0 else 'negative' if contrast['maximum']<0 else 'includes_zero'
    agreement=Counter(r['strict_numeric_contrast_sign_agreement'] for r in rows)
    screens={}
    for s in plan['screens']:
        sid=s['id'];bits=sum(r[sid+'_both_pass']<<i for i,r in enumerate(rows))
        screens[sid]=dict(pair_pass_bits=bits,passing_pairs=sum(r[sid+'_both_pass'] for r in rows),all_pairs_pass=bits==(1<<48)-1,any_pair_pass=bits!=0)
    return dict(triad_id=triad['triad_id'],models=triad['models'],mask=mask,sequence_method=method,mapping_definition=definition,
        sequence_robustness_key=[triad['triad_id'],mask,method],structural_robustness_key=[triad['triad_id'],mask,definition],
        pair_order_layout=[[p,f'{i:03b}'] for p in PERMUTATIONS for i in range(8)],pair_states=48,
        identical_correspondence_pairs=sum(r['identical_correspondence'] for r in rows),both_unique_numeric_pairs=unique,
        strict_numeric_contrast_sign_agreement_counts=dict(agreement),strict_all_pairs_contrast_difference_direction=direction,
        metric_difference_ranges=metrics,screens=screens)


def execute(plan,plan_path,source,bindings,out):
    config=out/'configuration.json';configuration=dict(plan_sha256=sha(plan_path),source_hashes=bindings)
    if config.exists():assert json.loads(config.read_text())==configuration,'Changed comparison checkpoint sources'
    else:
        with config.open('x') as f:f.write(json.dumps(configuration,indent=2)+'\n')
    checkpoint_root=out/'checkpoints';checkpoint_root.mkdir(exist_ok=True)
    checkpoint_bytes=sum(p.stat().st_size for p in checkpoint_root.glob('*.json.gz'))
    counts=Counter();seq_counts=Counter();joint_counts=Counter();passes=Counter();all_passes=Counter()
    maximum={f:0. for f in NUMERIC_FIELDS};groups=sgroups=states=seqrows=structrows=checked_reused=0;artifacts={};expected_checkpoints=set()
    names=['sequence_order_robustness.jsonl.gz','correspondence_pair_states.tsv.gz','correspondence_comparison_groups.jsonl.gz']
    with gzip.open(out/names[0],'wt',compresslevel=1) as sf,gzip.open(out/names[1],'wt',compresslevel=1) as pf,gzip.open(out/names[2],'wt',compresslevel=1) as jf:
        writer=csv.DictWriter(pf,fieldnames=pair_fields(plan),delimiter='\t',lineterminator='\n');writer.writeheader()
        for number,(triad,sequence,structural,st,tt) in enumerate(blocks(source),1):
            seqgroups=[];pairs=[];jointgroups=[]
            for mask in ['full','plddt70']:
                for method in METHODS:
                    rows=[sequence[mask,method,p] for p in PERMUTATIONS]
                    value=sequence_summary(triad,mask,method,rows,plan);seqgroups.append(value)
                    seq_counts[mask+':'+method+':'+value['strict_all_orders_contrast_direction']]+=1;sgroups+=1
                    for definition in DEFINITIONS:
                        combined=[]
                        for permutation in PERMUTATIONS:
                            for index in range(8):
                                s=sequence[mask,method,permutation];t=structural[mask,definition,index]
                                row=pair(triad,mask,method,permutation,definition,index,s,t,st[mask,method,permutation],tt[mask,definition,index],plan)
                                combined.append(row);states+=1
                                key=mask+':'+method+':'+definition;counts[key+':'+s['fit_status']+':'+t['fit_status']]+=1
                                for screen in plan['screens']:passes[key+':'+screen['id']]+=row[screen['id']+'_both_pass']
                                for field in NUMERIC_FIELDS:
                                    difference=row['sequence_minus_structural_'+field]
                                    if difference!='':maximum[field]=max(maximum[field],abs(difference))
                        value=joint_summary(triad,mask,method,definition,combined,plan);jointgroups.append(value);pairs.extend(combined);groups+=1
                        joint_counts[key+':'+value['strict_all_pairs_contrast_difference_direction']]+=1
                        for screen in plan['screens']:all_passes[key+':'+screen['id']]+=value['screens'][screen['id']]['all_pairs_pass']
            tid=triad['triad_id'];assert Path(tid).name==tid and tid not in ['', '.', '..']
            checkpoint=checkpoint_root/(tid+'.json.gz');expected_checkpoints.add(checkpoint.name)
            payload=dict(triad_id=tid,sequence_groups=seqgroups,pairs=pairs,comparison_groups=jointgroups)
            encoded=gzip.compress(json.dumps(payload,separators=(',',':')).encode(),compresslevel=1,mtime=0)
            if checkpoint.exists():assert checkpoint.read_bytes()==encoded,'False or changed comparison checkpoint';checked_reused+=1
            else:
                temporary=checkpoint.with_suffix('.tmp');temporary.write_bytes(encoded);temporary.replace(checkpoint)
                checkpoint_bytes+=len(encoded)
            artifacts[str(checkpoint.relative_to(out))]=sha(checkpoint)
            for value in seqgroups:sf.write(json.dumps(value,separators=(',',':'))+'\n')
            writer.writerows(pairs)
            for value in jointgroups:jf.write(json.dumps(value,separators=(',',':'))+'\n')
            seqrows+=len(sequence);structrows+=len(structural)
            if number%128==0:
                state=dict(status='running_full_triad_correspondence_comparison',completed_triads=number,total_triads=len(source['ready']),pair_states=states,checked_reused_triad_blocks=checked_reused)
                temp=out/'state.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(out/'state.json')
                assert checkpoint_bytes+sf.tell()+pf.tell()+jf.tell()<=plan['resources']['maximum_output_gib']*2**30
                assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
            if number%1000==0:print('Full correspondence comparison',number,'/',len(source['ready']),'pairs',states,flush=True)
    assert {p.name for p in checkpoint_root.glob('*.json.gz')}==expected_checkpoints
    summary=dict(ordered_model_triads=source['ordered_model_triads'],source_ready_triads=len(source['ready']),sequence_fit_rows=seqrows,structural_fit_rows=structrows,
        sequence_groups=sgroups,comparison_groups=groups,pair_states=states,pair_status_counts=dict(counts),sequence_direction_counts=dict(seq_counts),joint_direction_counts=dict(joint_counts),
        both_screen_pass_counts=dict(passes),all_pair_screen_pass_group_counts=dict(all_passes),maximum_absolute_metric_difference=maximum)
    assert all(summary[k]==v for k,v in plan['expected'].items())
    bind(bindings,config);verify(bindings)
    for name in names:artifacts[name]=sha(out/name)
    result=dict(status='complete_full_triad_correspondence_comparison_pending_independent_readback',**summary,plan_sha256=sha(plan_path),checked_reused_triad_blocks=checked_reused,
        source_hashes=bindings,artifacts=artifacts,scientific_eligibility=False,scope=plan['scope'])
    with (out/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','artifacts']}),flush=True);return result


def run(plan_path):
    plan=json.loads(Path(plan_path).read_text());assert shutil.disk_usage('.').free>=plan['resources']['minimum_free_disk_gib']*2**30
    source,bindings=load_sources(plan,plan_path);out=Path(plan['output']);out.mkdir(exist_ok=True,parents=True)
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (out/'receipt.json').exists(),'Completed output immutable'
        return execute(plan,plan_path,source,bindings,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);run(p.parse_args().plan)
