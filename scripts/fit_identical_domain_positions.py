#!/usr/bin/env python3
"""Fit identical complete domains by sequence offset, with full and joint-pLDDT masks."""
import argparse,csv,json,hashlib
from pathlib import Path
from functools import lru_cache
import numpy as np
from duplication_alignment_numeric_readback import load_pdb
from readback_cross_clan_alignments import rmsd
from assess_domain_alignment_geometry import geometry


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,rs):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args();p=json.loads(a.plan.read_text())
    for f,h in p['pins'].items():assert sha(f)==h
    out=Path(p['output']);assert not out.exists()
    candidates=[x for x in rows(p['candidates']) if x['all_complete_intervals_identical']=='1'];assert len(candidates)==48
    context={x['triad_key']:x for x in rows(p['context'])};pairs={};links=[]
    for c in candidates:
        used=set()
        for tk in json.loads(c['triad_keys_json']):
            t=context[tk];assert t['complete_intervals_identical']=='1';pair=(t['interval_a'],t['interval_b'])
            pk=hashlib.sha256(json.dumps(pair,separators=(',',':')).encode()).hexdigest();pairs[pk]=pair;used.add(pk)
        for pk in sorted(used):links.append({**{k:c[k] for k in ['family','gene_a','gene_b','pfam_accession','candidate_class','study_role','all_aligned_cores_identical']},'pair_key':pk})
    ids={x for pair in pairs.values() for x in pair};inputs={}
    for line in Path(p['inputs']).open():
        x=json.loads(line)
        if x['interval_id'] in ids and x['mask']=='full':inputs[x['interval_id']]=x
    assert set(inputs)==ids
    @lru_cache(maxsize=None)
    def data(iid):return load_pdb(inputs[iid])
    fits=[]
    for pk,(ia,ib) in sorted(pairs.items()):
        sa,xa,pa=data(ia);sb,xb,pb=data(ib);assert sa==sb;length=len(sa)
        for mask in ['full','joint_plddt70']:
            offsets=np.arange(length) if mask=='full' else np.flatnonzero((pa>=70)&(pb>=70));n=len(offsets)
            row=dict(pair_key=pk,interval_a=ia,interval_b=ib,mask=mask,interval_length=length,paired_residues=n,paired_fraction=n/length,sequence_offsets_zero_based=json.dumps(offsets.tolist(),separators=(',',':')),rmsd='',relative_rotation_curvature='',geometry_status='too_few_residues',mean_plddt_a='',mean_plddt_b='',n30_c70_pass=0)
            if n>=3:
                aa=xa[offsets];bb=xb[offsets];g=geometry(aa,bb);row.update(rmsd=rmsd(aa,bb),relative_rotation_curvature=g['relative_rotation_curvature'],geometry_status=g['geometry_status'],mean_plddt_a=float(pa[offsets].mean()),mean_plddt_b=float(pb[offsets].mean()),n30_c70_pass=int(n>=30 and 10*n>=7*length and g['geometry_status']=='unique_at_numeric_tolerance'))
            fits.append(row)
    out.mkdir(parents=True);save(out/'position_fits.tsv',fits);save(out/'candidate_pair_links.tsv',links)
    result=dict(status='complete_identical_domain_position_fits_pending_readback',plan_sha256=sha(a.plan),candidates=len(candidates),unique_interval_pairs=len(pairs),fits=len(fits),candidate_pair_links=len(links),masks=['full','joint_plddt70'],artifacts={x.name:sha(x) for x in out.iterdir()},scope='Exact sequence-offset pairing in identical complete domain intervals; full CA pairs and intersection of both pLDDT>=70 positions. These direct duplicate-pair RMSDs do not recalculate the three-protein reference asymmetry or establish biological structural changes. No structural alignment search, prediction or source substitution.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
