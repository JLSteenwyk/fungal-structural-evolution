#!/usr/bin/env python3
"""Independently verify sequence-offset fits via quaternion rotations."""
import argparse,csv,json,hashlib,math,statistics
from pathlib import Path
from functools import lru_cache
import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from readback_domain_triad_common_fits import quaternion_fit


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();p=json.loads(a.plan.read_text());root=Path(p['output']);r=json.loads((root/'receipt.json').read_text());assert r['plan_sha256']==sha(a.plan)
    for f,h in p['pins'].items():assert sha(f)==h
    for f,h in r['artifacts'].items():assert sha(root/f)==h
    candidates=[x for x in read(p['candidates']) if x['all_complete_intervals_identical']=='1'];context={x['triad_key']:x for x in read(p['context'])};pairs={};expectedlinks=[]
    for c in candidates:
        local=set()
        for tk in json.loads(c['triad_keys_json']):
            t=context[tk];pair=[t['interval_a'],t['interval_b']];pk=hashlib.sha256(json.dumps(pair,separators=(',',':')).encode()).hexdigest();pairs[pk]=pair;local.add(pk)
        for pk in local:expectedlinks.append(tuple([c[k] for k in ['family','gene_a','gene_b','pfam_accession','candidate_class','study_role','all_aligned_cores_identical']]+[pk]))
    links=read(root/'candidate_pair_links.tsv');fields=['family','gene_a','gene_b','pfam_accession','candidate_class','study_role','all_aligned_cores_identical','pair_key'];assert len(links)==len(expectedlinks)==r['candidate_pair_links'] and {tuple(x[k] for k in fields) for x in links}==set(expectedlinks)
    ids={i for pair in pairs.values() for i in pair};inputs={}
    for line in Path(p['inputs']).open():
        x=json.loads(line)
        if x['interval_id'] in ids and x['mask']=='full':inputs[x['interval_id']]=x
    @lru_cache(maxsize=None)
    def coords(iid):
        x=inputs[iid];assert sha(x['path'])==x['sha256'];letters=[];positions=[];xyz=[];confidence=[]
        for line in Path(x['path']).read_text().splitlines():
            if line.startswith('ATOM  '):
                assert line[12:16].strip()=='CA';letters.append(protein_letters_3to1[line[17:20].strip().title()]);positions.append(int(line[22:26]));xyz.append([float(line[30:38]),float(line[38:46]),float(line[46:54])]);confidence.append(float(line[60:66]))
        assert ''.join(letters)==x['sequence'] and positions==x['original_positions'];return ''.join(letters),np.array(xyz),confidence
    fits=read(root/'position_fits.tsv');seen=set();maximum=0.;summaries={}
    for x in fits:
        key=x['pair_key'],x['mask'];assert key not in seen;seen.add(key);ia,ib=pairs[key[0]];assert x['interval_a']==ia and x['interval_b']==ib
        sa,xa,pa=coords(ia);sb,xb,pb=coords(ib);assert sa==sb
        offsets=[i for i in range(len(sa)) if x['mask']=='full' or (pa[i]>=70 and pb[i]>=70)];n=len(offsets)
        assert json.loads(x['sequence_offsets_zero_based'])==offsets and int(x['paired_residues'])==n and int(x['interval_length'])==len(sa)
        assert math.isclose(float(x['paired_fraction']),n/len(sa),abs_tol=1e-12)
        if n<3:
            assert x['geometry_status']=='too_few_residues' and all(x[f]=='' for f in ['rmsd','relative_rotation_curvature','mean_plddt_a','mean_plddt_b']) and x['n30_c70_pass']=='0';continue
        distance,curvature,status=quaternion_fit(xa[offsets],xb[offsets])
        assert x['geometry_status']==status
        for name,wanted in [('rmsd',distance),('relative_rotation_curvature',curvature),('mean_plddt_a',sum(pa[i] for i in offsets)/n),('mean_plddt_b',sum(pb[i] for i in offsets)/n)]:assert math.isclose(float(x[name]),wanted,abs_tol=1e-9,rel_tol=1e-9),(key,name)
        assert int(x['n30_c70_pass'])==int(n>=30 and n/len(sa)>=.7 and status=='unique_at_numeric_tolerance')
        maximum=max(maximum,abs(distance-float(x['rmsd'])))
    assert seen=={(pk,m) for pk in pairs for m in ['full','joint_plddt70']} and len(fits)==r['fits'] and len(pairs)==r['unique_interval_pairs'] and len(candidates)==r['candidates']
    for mask in ['full','joint_plddt70']:
        xs=[x for x in fits if x['mask']==mask];values=[float(x['rmsd']) for x in xs if x['rmsd']!=''];summaries[mask]=dict(fits=len(xs),computed=len(values),n30_c70_passing=sum(x['n30_c70_pass']=='1' for x in xs),rmsd_min=min(values),rmsd_median=statistics.median(values),rmsd_max=max(values))
    result=dict(status='passed_full_identical_domain_position_fit_readback',candidates=len(candidates),interval_pairs=len(pairs),fits=len(fits),pdbs_checked=len(ids),maximum_rmsd_difference=maximum,mask_summaries=summaries,plan_sha256=sha(a.plan),producer_receipt_sha256=sha(root/'receipt.json'),checker_sha256=sha(__file__),quaternion_helper_sha256=sha(Path(__file__).with_name('readback_domain_triad_common_fits.py')),scope='Exact pair/link universe, every sequence offset/mask, raw PDB hash/sequence/coordinates, RMSD/geometry/confidence and screen decisions checked. Quaternion rotation independent of producer SVD. Pair distances do not verify three-protein directional asymmetry or prediction accuracy.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
