#!/usr/bin/env python3
"""Constrain identical duplicate domains to sequence offsets in shared reference cores."""
import argparse,csv,gzip,hashlib,itertools,json
from pathlib import Path
from functools import lru_cache
from collections import Counter,defaultdict
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from fit_domain_triad_common_residues import fit_triplet


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
    for f,h in p['pins'].items():assert sha(f)==h
    out=Path(p['output']);assert not out.exists()
    candidates=[x for x in read(p['candidates']) if x['all_complete_intervals_identical']=='1'];assert len(candidates)==48
    needed={t for c in candidates for t in json.loads(c['triad_keys_json'])};triads={}
    for x in map(json.loads,Path(p['triads']).open()):
        if x['triad_key'] in needed:triads[x['triad_key']]=[x['interval_'+role] for role in ['a','b','reference']]
    assert set(triads)==needed
    ids={i for ts in triads.values() for i in ts};inputs={}
    for x in map(json.loads,Path(p['inputs']).open()):
        if x['interval_id'] in ids and x['mask']=='full':inputs[x['interval_id']]=x
    assert set(inputs)==ids
    @lru_cache(maxsize=None)
    def data(iid):
        r=inputs[iid];s,x,p=load_pdb(r);assert r['original_positions']==list(range(r['start'],r['end']+1));return s,x,p
    sample=fit_triplet([np.eye(3)]*3,['AAA']*3,[np.full(3,90.)]*3);fits=[];seen=set();bytriad=defaultdict(list)
    with gzip.open(p['maps'],'rt') as f:
        for line in f:
            m=json.loads(line);tk=m['triad_key']
            if tk not in needed:continue
            assert not any(m['edge_exclusions'])
            ids3=triads[tk];loaded=[data(i) for i in ids3];starts=[inputs[i]['start'] for i in ids3];lengths=[len(x[0]) for x in loaded];assert lengths==m['original_interval_lengths'] and loaded[0][0]==loaded[1][0]
            for definition,anchor in itertools.product(['reference_common','cycle_consistent'],['a','b']):
                source=m[definition+'_triples'];chosen=0 if anchor=='a' else 1;triples=[]
                for t in source:
                    offset=t[chosen]-starts[chosen];rpos=t[2]-starts[2];assert 0<=offset<lengths[0] and 0<=rpos<lengths[2]
                    if m['mask']=='plddt70' and not all(loaded[k][2][j]>=70 for k,j in enumerate([offset,offset,rpos])):continue
                    triples.append([starts[0]+offset,starts[1]+offset,t[2]])
                key=(tk,m['mask'],*m['orders'],definition,anchor);assert key not in seen;seen.add(key)
                n=len(triples);row=dict(triad_key=tk,mask=m['mask'],order_ab=m['orders'][0],order_ar=m['orders'][1],order_br=m['orders'][2],mapping_definition=definition,reference_anchor=anchor,source_common_residues=len(source),common_residues=n,removed_for_joint_confidence=len(source)-n,control_triples_json=json.dumps(triples,separators=(',',':')))
                for k,role in enumerate(['a','b','reference']):row['length_'+role]=lengths[k];row['coverage_'+role]=n/lengths[k]
                row.update({field:'' for field in sample});row['fit_status']='fewer_than_three_common_residues';row['n30_c70_pass']=0
                if n>=3:
                    indices=[[t[k]-starts[k] for t in triples] for k in range(3)];coords=[loaded[k][1][indices[k]] for k in range(3)];letters=[''.join(loaded[k][0][j] for j in indices[k]) for k in range(3)];conf=[loaded[k][2][indices[k]] for k in range(3)]
                    row.update(fit_triplet(coords,letters,conf));assert row['sequence_identity_ab']==1
                    row['n30_c70_pass']=int(n>=30 and all(10*n>=7*l for l in lengths) and row['fit_status']=='computed_unique_at_numeric_tolerance')
                fits.append(row);bytriad[tk].append(row)
    assert set(bytriad)==needed and all(len(xs)==64 for xs in bytriad.values())
    summaries=[]
    for c in candidates:
        pool=[x for t in json.loads(c['triad_keys_json']) for x in bytriad[t]];eligible=[x for x in pool if x['n30_c70_pass']];complete=len(eligible)==len(pool)
        low=min((x['rmsd_ar_minus_br'] for x in eligible),default=None);high=max((x['rmsd_ar_minus_br'] for x in eligible),default=None)
        if not complete:label='incomplete'
        elif low>.1:label='a_farther_from_reference'
        elif high<-.1:label='b_farther_from_reference'
        elif low<-.1 and high>.1:label='direction_flip_beyond_margin'
        else:label='touches_or_enters_margin_band'
        summaries.append(dict(c,control_expected_fits=len(pool),control_passing_fits=len(eligible),control_complete=int(complete),control_contrast_min='' if low is None else low,control_contrast_max='' if high is None else high,control_direction_margin_0_1=label,control_retains_original_direction=int(label==c['mafft_structural_direction'])))
    out.mkdir(parents=True);save(out/'sequence_locked_fits.tsv',fits);save(out/'candidate_control_summary.tsv',summaries)
    result=dict(status='complete_sequence_locked_reference_cores_pending_readback',plan_sha256=sha(a.plan),candidates=len(candidates),triads=len(needed),fits=len(fits),passing_fits=sum(x['n30_c70_pass'] for x in fits),candidate_directions=dict(Counter(x['control_direction_margin_0_1'] for x in summaries)),retains_original_direction=sum(x['control_retains_original_direction'] for x in summaries),artifacts={f.name:sha(f) for f in out.iterdir()},scope='Original shared reference positions retained; use each duplicate reference mapping separately, force identical duplicate sequence offsets, then apply joint three-protein confidence for pLDDT70. Both mapping definitions, eight orders, two masks and two anchors retained (64 per triad). Conditional correspondence sensitivity, not ancestral change or prediction-error calibration.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
