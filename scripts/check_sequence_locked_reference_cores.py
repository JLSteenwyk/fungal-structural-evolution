#!/usr/bin/env python3
"""Independently reconstruct all constrained reference cores and quaternion fits."""
import argparse,csv,gzip,json,hashlib,itertools,math
from collections import Counter,defaultdict
from pathlib import Path
from functools import lru_cache
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from readback_domain_triad_common_fits import core_metrics


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();p=json.loads(a.plan.read_text());root=Path(p['output']);receipt=json.loads((root/'receipt.json').read_text());assert receipt['plan_sha256']==sha(a.plan)
    for f,h in p['pins'].items():assert sha(f)==h
    for f,h in receipt['artifacts'].items():assert sha(root/f)==h
    candidates=[x for x in read(p['candidates']) if x['all_complete_intervals_identical']=='1'];wanted={t for c in candidates for t in json.loads(c['triad_keys_json'])}
    triads={r['triad_key']:[r['interval_'+role] for role in ['a','b','reference']] for r in map(json.loads,Path(p['triads']).open()) if r['triad_key'] in wanted};ids=set(itertools.chain.from_iterable(triads.values()));inputs={}
    for r in map(json.loads,Path(p['inputs']).open()):
        if r['interval_id'] in ids and r['mask']=='full':inputs[r['interval_id']]=r
    maps={}
    with gzip.open(p['maps'],'rt') as f:
        for r in map(json.loads,f):
            if r['triad_key'] in wanted:
                k=(r['triad_key'],r['mask'],*map(str,r['orders']));assert k not in maps;maps[k]=r
    @lru_cache(maxsize=None)
    def data(i):return load_pdb(inputs[i])
    fits=read(root/'sequence_locked_fits.tsv');seen=set();bytriad=defaultdict(list);maxerr=0.
    for r in fits:
        key=tuple(r[k] for k in ['triad_key','mask','order_ab','order_ar','order_br','mapping_definition','reference_anchor']);assert key not in seen;seen.add(key)
        src=maps[key[:5]];assert not any(src['edge_exclusions']);raw=src[r['mapping_definition']+'_triples'];iids=triads[r['triad_key']];loaded=[data(i) for i in iids];starts=np.array([inputs[i]['start'] for i in iids]);lengths=[len(x[0]) for x in loaded];assert loaded[0][0]==loaded[1][0]
        relative=np.array(raw,dtype=int).reshape(-1,3)-starts
        if r['reference_anchor']=='a':relative[:,1]=relative[:,0]
        else:assert r['reference_anchor']=='b';relative[:,0]=relative[:,1]
        assert all(((relative[:,k]>=0)&(relative[:,k]<lengths[k])).all() for k in range(3))
        if r['mask']=='plddt70':relative=relative[np.logical_and.reduce([loaded[k][2][relative[:,k]]>=70 for k in range(3)])]
        else:assert r['mask']=='full'
        assert json.loads(r['control_triples_json'])==(relative+starts).tolist();n=len(relative)
        assert int(r['source_common_residues'])==len(raw) and int(r['common_residues'])==n and int(r['removed_for_joint_confidence'])==len(raw)-n
        for k,role in enumerate(['a','b','reference']):assert int(r['length_'+role])==lengths[k] and math.isclose(float(r['coverage_'+role]),n/lengths[k],abs_tol=1e-12)
        if n>=3:
            coords=[loaded[k][1][relative[:,k]] for k in range(3)];letters=[''.join(loaded[k][0][j] for j in relative[:,k]) for k in range(3)];conf=[loaded[k][2][relative[:,k]] for k in range(3)];metrics=core_metrics(coords,letters,conf);assert metrics['sequence_identity_ab']==1
            for field,v in metrics.items():
                if isinstance(v,str):assert r[field]==v,(key,field)
                else:
                    assert math.isclose(float(r[field]),v,rel_tol=1e-9,abs_tol=1e-9),(key,field)
                    if field.startswith('rmsd'):maxerr=max(maxerr,abs(float(r[field])-v))
            passing=int(n>=30 and min(n/l for l in lengths)>=.7 and metrics['fit_status']=='computed_unique_at_numeric_tolerance')
        else:
            assert r['fit_status']=='fewer_than_three_common_residues';passing=0
            for f in ['rmsd_ab','rmsd_ar','rmsd_br','rmsd_ar_minus_br']:assert r[f]==''
        assert int(r['n30_c70_pass'])==passing;bytriad[r['triad_key']].append(r)
    expected={(tk,mask,*map(str,orders),definition,anchor) for tk,mask,orders,definition,anchor in itertools.product(wanted,['full','plddt70'],list(itertools.product([0,1],repeat=3)),['reference_common','cycle_consistent'],['a','b'])};assert seen==expected
    keys=['family','gene_a','gene_b','pfam_accession'];summary=read(root/'candidate_control_summary.tsv');index={tuple(x[k] for k in keys):x for x in summary};assert len(index)==len(summary)==len(candidates)==receipt['candidates'];counts=Counter();retained=0;changed=[]
    for c in candidates:
        s=index[tuple(c[k] for k in keys)]
        for f,v in c.items():assert s[f]==v
        pool=[r for t in json.loads(c['triad_keys_json']) for r in bytriad[t]];passing=[r for r in pool if r['n30_c70_pass']=='1'];complete=len(pool)==len(passing)
        assert int(s['control_expected_fits'])==len(pool) and int(s['control_passing_fits'])==len(passing) and int(s['control_complete'])==int(complete)
        low=min((float(r['rmsd_ar_minus_br']) for r in passing),default=None);high=max((float(r['rmsd_ar_minus_br']) for r in passing),default=None)
        for f,v in [('control_contrast_min',low),('control_contrast_max',high)]:assert (s[f]=='' if v is None else math.isclose(float(s[f]),v,abs_tol=1e-12,rel_tol=1e-12))
        label='incomplete' if not complete else 'a_farther_from_reference' if low>.1 else 'b_farther_from_reference' if high<-.1 else 'direction_flip_beyond_margin' if low<-.1 and high>.1 else 'touches_or_enters_margin_band'
        assert s['control_direction_margin_0_1']==label and int(s['control_retains_original_direction'])==int(label==c['mafft_structural_direction']);counts[label]+=1;retained+=int(label==c['mafft_structural_direction'])
        if label!=c['mafft_structural_direction']:changed.append({k:s[k] for k in keys+['pfam_name','candidate_class','all_aligned_cores_identical','control_direction_margin_0_1','control_contrast_min','control_contrast_max']})
    assert dict(counts)==receipt['candidate_directions'] and retained==receipt['retains_original_direction'] and len(fits)==receipt['fits'] and sum(x['n30_c70_pass']=='1' for x in fits)==receipt['passing_fits']
    result=dict(status='passed_full_sequence_locked_reference_core_readback',candidates=len(candidates),triads=len(wanted),fits=len(fits),maximum_rmsd_difference=maxerr,candidate_directions=dict(counts),retains_original_direction=retained,changed_candidates=changed,plan_sha256=sha(a.plan),producer_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),quaternion_helper_sha256=sha(Path(__file__).with_name('readback_domain_triad_common_fits.py')),scope='All reference triples independently reconstructed by relative-coordinate arrays; every mask/count/metric checked with quaternion fits, full alternative grid and candidate classifications verified. Conditional correspondence sensitivity only; no prediction accuracy or evolutionary causality claim.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
