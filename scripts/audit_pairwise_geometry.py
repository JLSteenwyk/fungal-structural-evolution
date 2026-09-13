#!/usr/bin/env python3
"""Read back every pair/threshold eligibility and independently spot-check geometry."""
import argparse,csv,gzip,itertools,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import pdist
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--comparisons',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new audit output')
    checked_receipt(a.mapping);checked_receipt(a.comparisons)
    receipt=json.loads((a.comparisons/'receipt.json').read_text())
    if receipt['mapping_receipt_sha256']!=sha(a.mapping/'receipt.json'):raise ValueError('Mapping identity differs')
    mr=json.loads((a.mapping/'receipt.json').read_text());matrix=ROOT/'results/phylogeny/profile-matrix-50-v1'
    checked_receipt(matrix)
    if sha(matrix/'receipt.json')!=mr['matrix_receipt_sha256']:raise ValueError('Matrix identity differs')
    aa={r.id:np.array(list(str(r.seq))) for r in SeqIO.parse(matrix/'matrix.faa','fasta')}
    links=read_table(a.mapping/'marker_structure_links.tsv');by={};groups=defaultdict(list);mapped=defaultdict(dict)
    for r in links:
        key=(r['marker'],r['taxon_id'])
        if key in by:raise ValueError('Duplicate marker/taxon mapping')
        by[key]=r;groups[key[0]].append(key[1])
    with gzip.open(a.mapping/'matrix_to_structure_residues.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=(r['marker'],r['taxon_id']);col=int(r['matrix_column_1based'])-1
            if col in mapped[key]:raise ValueError('Duplicate mapped position')
            mapped[key][col]=(int(r['protein_residue_1based']),float(r['ca_plddt']))
    expected={(m,x,y) for m,ts in groups.items() for x,y in itertools.combinations(sorted(ts),2)}
    seen={};counts=Counter();samples={};cache_key=None;cache=None
    for name,accepted in [('pairwise_metrics.tsv',True),('pairwise_exclusions.tsv',False)]:
        with (a.comparisons/name).open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                key=(r['marker'],r['taxon_a'],r['taxon_b']);cut=int(r['plddt_cutoff'])
                if key not in expected or cut not in (50,70,90):raise ValueError('Unexpected pair/threshold')
                bit={50:1,70:2,90:4}[cut]
                if seen.get(key,0)&bit:raise ValueError('Duplicate pair/threshold')
                seen[key]=seen.get(key,0)|bit
                if key!=cache_key:
                    m,x,y=key;ma,mb=mapped[m,x],mapped[m,y];cols=np.array(sorted(ma.keys()&mb.keys()),dtype=int)
                    confidence=np.array([min(ma[c][1],mb[c][1]) for c in cols]);canonical=np.isin(aa[x][cols],list('ACDEFGHIKLMNPQRSTVWY'))&np.isin(aa[y][cols],list('ACDEFGHIKLMNPQRSTVWY'))
                    cache_key=key;cache=(cols,confidence,canonical)
                cols,confidence,canonical=cache;selected=cols[(confidence>=cut)&canonical];n=len(selected);shared=len(cols)
                ok=n>=50 and 2*n>=shared
                if (ok!=accepted or n!=int(r['compared_residues']) or shared!=int(r['shared_profile_positions']) or not math.isclose(float(r['fraction_of_shared_positions_compared']),n/shared if shared else 0,abs_tol=1e-12)):raise ValueError('Eligibility readback failed')
                m,x,y=key
                if (r['same_model_coordinates']=='True')!=(by[m,x]['model_sha256']==by[m,y]['model_sha256']):raise ValueError('Model identity differs')
                counts[f'{cut}_{"accepted" if accepted else "excluded"}']+=1
                if accepted:
                    diff=float(np.mean(aa[x][selected]!=aa[y][selected]))
                    if not math.isclose(diff,float(r['uncorrected_sequence_difference']),abs_tol=1e-12):raise ValueError('Sequence difference differs')
                    for field in ['ca_superposition_rmsd_angstrom','local_distance_mean_absolute_change_angstrom','local_distance_rms_change_angstrom']:
                        if r[field] and (not math.isfinite(float(r[field])) or float(r[field])<0):raise ValueError('Invalid geometry metric')
                    rank=sha_text('|'.join(map(str,(*key,cut))))
                    # Deterministic smallest hash per marker and cutoff avoids choosing by outcome.
                    sk=(m,cut)
                    if sk not in samples or rank<samples[sk][0]:samples[sk]=(rank,dict(r),selected.copy())
    if set(seen)!=expected or any(v!=7 for v in seen.values()):raise ValueError('Incomplete pair/threshold grid')
    if sum(v for k,v in counts.items() if k.endswith('_accepted'))!=receipt['comparison_rows'] or sum(v for k,v in counts.items() if k.endswith('_excluded'))!=receipt['excluded_rows']:raise ValueError('Receipt counts differ')
    checked=[];models={}
    for _,r,cols in samples.values():
        m,x,y=r['marker'],r['taxon_a'],r['taxon_b'];positions=[];coords=[]
        for t in (x,y):
            link=by[m,t];path=ROOT/link['model_path']
            if str(path) not in models:
                if sha(path)!=link['model_sha256']:raise ValueError('Changed model')
                c=MMCIF2Dict(str(path));models[str(path)]={int(p):[float(x),float(y),float(z)] for atom,p,x,y,z in zip(c['_atom_site.label_atom_id'],c['_atom_site.label_seq_id'],c['_atom_site.Cartn_x'],c['_atom_site.Cartn_y'],c['_atom_site.Cartn_z']) if atom=='CA'}
            pos=np.array([mapped[m,t][int(c)][0] for c in cols]);positions.append(pos);coords.append(np.array([models[str(path)][p] for p in pos]))
        u,v=coords;uc=u-u.mean(0);vc=v-v.mean(0);rot,_=Rotation.align_vectors(vc,uc)
        rms=np.sqrt(np.mean(np.sum((rot.apply(uc)-vc)**2,axis=1)))
        dx,dy=pdist(u),pdist(v);px,py=positions;ii,jj=np.triu_indices(len(u),1)
        keep=((dx<=15)|(dy<=15))&(np.abs(px[ii]-px[jj])>=3)&(np.abs(py[ii]-py[jj])>=3)
        delta=np.abs(dx[keep]-dy[keep]);values={'ca_superposition_rmsd_angstrom':rms,'local_distance_pairs':len(delta)}
        if len(delta):values.update(local_distance_mean_absolute_change_angstrom=float(delta.mean()),local_distance_rms_change_angstrom=float(np.sqrt(np.mean(delta**2))))
        for field,val in values.items():
            if not math.isclose(val,float(r[field]),rel_tol=1e-8,abs_tol=1e-8):raise ValueError('Independent geometry mismatch: '+field)
        checked.append({k:r[k] for k in ['marker','taxon_a','taxon_b','plddt_cutoff']})
    a.output.mkdir(parents=True)
    result={'status':'complete_pairwise_grid_and_sampled_geometry_readback','mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'comparison_receipt_sha256':sha(a.comparisons/'receipt.json'),'script_sha256':sha(Path(__file__)),'expected_pairs':len(expected),'all_pair_threshold_rows':sum(counts.values()),'counts':dict(counts),'independent_geometry_checks':len(checked),'geometry_selection':'Smallest SHA256 pair/threshold identity per accepted marker and threshold; SciPy Rotation and condensed pairwise distances, no production geometry helper.','geometry_checked_rows':checked,'interpretation':'Every eligibility, compared-site count and accepted sequence difference verified; geometry numerically checked on a deterministic subset only. Focal confidence only; PAE and six-residue context remain separate gates. No evolutionary-rate validation claimed.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='geometry_checked_rows'},indent=2))


def sha_text(s):
    import hashlib
    return hashlib.sha256(s.encode()).hexdigest()


if __name__=='__main__':main()
