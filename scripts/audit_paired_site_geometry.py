#!/usr/bin/env python3
"""Read back complete paired-site grids and independently recompute sampled geometry."""
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['comparisons','inputs','snapshot','pae','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable audit output')
    receipts={n:checked_receipt(getattr(a,n)) for n in ['comparisons','inputs','snapshot','pae']}
    for n in ['inputs','snapshot','pae']:
        if receipts['comparisons']['source_receipts'][n]['sha256']!=sha(getattr(a,n)/'receipt.json'):raise ValueError('Source lineage mismatch')
    ready={r['marker']:r for r in read_table(a.inputs/'marker_summary.tsv') if r['status']=='ready_for_inference'}
    aa={};states={};expected=set();samples={};seen=set();counts=Counter()
    for m in ready:
        aa[m]={r.id:np.array(list(str(r.seq))) for r in SeqIO.parse(a.inputs/m/'aa.faa','fasta')}
        states[m]={r.id:np.array(list(str(r.seq))) for r in SeqIO.parse(a.inputs/m/'3di.faa','fasta')}
        if set(aa[m])!=set(states[m]):raise ValueError('Paired taxon sets differ')
        for t in aa[m]:
            if not np.array_equal(aa[m][t]=='?',states[m][t]=='?'):raise ValueError('Paired masks differ')
        expected.update((m,x,y) for x,y in combinations(sorted(aa[m]),2))
    for file,accepted in [('paired_site_geometry.tsv',True),('coverage_exclusions.tsv',False)]:
        for r in read_table(a.comparisons/file):
            key=(r['marker'],r['taxon_a'],r['taxon_b']);m,x,y=key
            if key not in expected or key in seen:raise ValueError('Unexpected/duplicate pair')
            seen.add(key);u,v=aa[m][x],aa[m][y];common=(u!='?')&(v!='?');n=int(common.sum());na=int((u!='?').sum());nb=int((v!='?').sum())
            if accepted!=(n>=50 and 2*n>=max(na,nb)):raise ValueError('Eligibility mismatch')
            for field,value in [('tree_taxa',len(aa[m])),('tree_alignment_columns',len(u)),('observed_a',na),('observed_b',nb),('geometry_matched_residues',n)]:
                if int(r[field])!=value:raise ValueError('Dimension mismatch '+field)
            if not math.isclose(float(r['fraction_of_tree_columns']),n/len(u),abs_tol=1e-12):raise ValueError('Coverage fraction mismatch')
            counts['accepted' if accepted else 'excluded']+=1
            if accepted:
                for field,value in [('uncorrected_aa_difference',np.mean(u[common]!=v[common])),('uncorrected_3di_difference',np.mean(states[m][x][common]!=states[m][y][common]))]:
                    if not math.isclose(float(r[field]),value,abs_tol=1e-12):raise ValueError('Character difference mismatch')
                rank=hashlib.sha256('|'.join(key).encode()).hexdigest()
                if m not in samples or rank<samples[m][0]:samples[m]=(rank,r,np.flatnonzero(common))
    if seen!=expected or counts['accepted']!=receipts['comparisons']['accepted_taxon_pairs'] or counts['excluded']!=receipts['comparisons']['excluded_taxon_pairs']:raise ValueError('Grid/receipt count mismatch')
    links={(r['marker'],r['taxon_id']):r for r in read_table(a.snapshot/'marker_structure_links.tsv')}
    mapped=defaultdict(dict)
    with gzip.open(a.snapshot/'matrix_to_structure_residues.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):mapped[r['marker'],r['taxon_id']][int(r['matrix_column_1based'])]=int(r['protein_residue_1based'])
    pae={(r['model_id'],str(r['version'])):r for r in json.loads((a.pae/'pae_manifest.json').read_text())};checks=[]
    for m,(_,r,common) in sorted(samples.items()):
        cols=[int(r['matrix_column_1based']) for r in read_table(a.inputs/m/'columns.tsv')];coords=[];positions=[];errors=[]
        for t in (r['taxon_a'],r['taxon_b']):
            link=links[m,t];path=ROOT/link['model_path']
            if sha(path)!=link['model_sha256']:raise ValueError('Changed coordinates')
            cif=MMCIF2Dict(str(path));ca={}
            for atom,pos,x,y,z in zip(cif['_atom_site.label_atom_id'],cif['_atom_site.label_seq_id'],cif['_atom_site.Cartn_x'],cif['_atom_site.Cartn_y'],cif['_atom_site.Cartn_z']):
                if atom=='CA':
                    if int(pos) in ca:raise ValueError('Duplicate CA')
                    ca[int(pos)]=[float(x),float(y),float(z)]
            pos=np.array([mapped[m,t][cols[i]] for i in common]);positions.append(pos);coords.append(np.array([ca[i] for i in pos]))
            pr=pae[link['model_id'],link['model_version']];compressed=(ROOT/pr['path']).read_bytes();raw=gzip.decompress(compressed)
            if hashlib.sha256(compressed).hexdigest()!=pr['gzip_sha256'] or hashlib.sha256(raw).hexdigest()!=pr['json_sha256']:raise ValueError('PAE checksum mismatch')
            data=json.loads(raw);data=data[0] if isinstance(data,list) else data;mat=np.array(data['predicted_aligned_error'],dtype=float)
            if mat.shape!=(int(pr['length']),int(pr['length'])) or not np.isfinite(mat).all() or (mat<0).any():raise ValueError('Invalid PAE matrix')
            errors.append(mat)
        x,y=coords;xc=x-x.mean(0);yc=y-y.mean(0);rot,_=Rotation.align_vectors(yc,xc)
        rms=float(np.sqrt(np.mean(np.sum((rot.apply(xc)-yc)**2,axis=1))))
        dx,dy=pdist(x),pdist(y);i,j=np.triu_indices(len(x),1);px,py=positions
        local=((dx<=15)|(dy<=15))&(abs(px[i]-px[j])>=3)&(abs(py[i]-py[j])>=3)
        confident=local.copy()
        for mat,pos in zip(errors,positions):confident&=(mat[pos[i]-1,pos[j]-1]<=10)&(mat[pos[j]-1,pos[i]-1]<=10)
        delta=abs(dx-dy)
        values={'ca_superposition_rmsd_angstrom':rms,'local_distance_pairs':int(local.sum()),'pae10_local_pairs':int(confident.sum())}
        if local.any():values.update(local_distance_mean_absolute_change_angstrom=float(delta[local].mean()),local_distance_rms_change_angstrom=float(np.sqrt(np.mean(delta[local]**2))),pae10_local_retained_fraction=float(confident.sum()/local.sum()))
        if confident.any():values['pae10_local_mean_absolute_change_angstrom']=float(delta[confident].mean())
        for field,value in values.items():
            if not math.isclose(float(r[field]),value,rel_tol=1e-8,abs_tol=1e-8):raise ValueError('Independent geometry mismatch '+field)
        checks.append({k:r[k] for k in ['marker','taxon_a','taxon_b']})
    a.output.mkdir(parents=True)
    result={'status':'passed_complete_paired_grid_character_and_sampled_geometry_readback','source_receipts':{n:sha(getattr(a,n)/'receipt.json') for n in receipts},'script_sha256':sha(Path(__file__)),'pair_rows':len(seen),'counts':dict(counts),'independent_geometry_checks':len(checks),'sample_selection':'Smallest SHA256 of marker and taxon pair per accepted marker, independent of metric values','checked_pairs':checks,'interpretation':'All pair identities, shared paired masks, eligibility, dimensions and AA/3Di differences checked. Independent SciPy rotation, condensed distances and direct directional PAE indexing checked on one pair per accepted marker only. Not a full numerical audit of all geometry or a revalidation of upstream confidence masks, biological orthology, branch additivity or rate inference.'}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='checked_pairs'},indent=2))

if __name__=='__main__':main()
