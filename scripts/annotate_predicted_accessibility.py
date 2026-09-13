#!/usr/bin/env python3
"""Compute residue solvent accessibility on frozen isolated predicted protein chains."""
import argparse,csv,fcntl,gzip,hashlib,inspect,json,math,os,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import Bio,numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley,ATOMIC_RADII
from Bio.Data.IUPACData import protein_letters_3to1
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha
THREE={k.upper():v for k,v in protein_letters_3to1.items()}


def run_one(job):
    model,output,config_hash,points=job;out=Path(output);sid=model['model_id'];rp=out/(sid+'.receipt.json');table=out/(sid+'.residues.tsv.gz')
    if rp.exists():
        r=json.loads(rp.read_text())
        if r['config_sha256']!=config_hash or r['table_sha256']!=sha(table):raise ValueError('Changed accessibility output')
        return r
    start=time.monotonic();path=ROOT/model['path']
    if sha(path)!=model['sha256']:raise ValueError('Changed predicted coordinates')
    structure=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False).get_structure(sid,path)
    ms=list(structure.get_models())
    if len(ms)!=1 or len(list(ms[0].get_chains()))!=1:raise ValueError('Expected one isolated predicted chain')
    residues=list(ms[0].get_residues());seq=''.join(THREE.get(r.resname,'?') for r in residues)
    if len(seq)!=model['length'] or hashlib.sha256(seq.encode()).hexdigest()!=model['sequence_sha256']:raise ValueError('Residue sequence differs')
    if [r.id[1] for r in residues]!=list(range(1,len(seq)+1)):raise ValueError('Residue index grid differs')
    atoms=list(ms[0].get_atoms())
    if any(a.is_disordered() or a.element not in {'C','N','O','S'} or not np.isfinite(a.coord).all() for a in atoms):raise ValueError('Unexpected/disordered/nonfinite predicted atom')
    if any(r.id[0]!=' ' or r.id[2]!=' ' or not r.has_id('CA') for r in residues):raise ValueError('Nonstandard or incomplete residue identity')
    sr=ShrakeRupley(probe_radius=1.4,n_points=points);sr.compute(ms[0],level='R')
    rows=[]
    for i,r in enumerate(residues,1):
        asa=float(r.sasa);plddt=float(r['CA'].bfactor)
        if not math.isfinite(asa) or asa<0 or not math.isfinite(plddt) or not 0<=plddt<=100:raise ValueError('Invalid accessibility/confidence')
        if not math.isclose(asa,sum(float(a.sasa) for a in r),rel_tol=1e-10,abs_tol=1e-8):raise ValueError('Residue/atom ASA sum differs')
        rows.append({'model_id':sid,'sequence_sha256':model['sequence_sha256'],'protein_residue_1based':i,'amino_acid':seq[i-1],'sasa_angstrom_squared':asa,'ca_plddt':plddt,'heavy_atom_count':len(r),'context':'isolated_predicted_chain'})
    partial=table.with_suffix('.partial')
    with gzip.open(partial,'wt') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    partial.replace(table)
    result={'status':'complete_predicted_chain_accessibility','model_id':sid,'model_sha256':model['sha256'],'sequence_sha256':model['sequence_sha256'],'config_sha256':config_hash,'residues':len(rows),'heavy_atoms':len(atoms),'total_sasa_angstrom_squared':sum(r['sasa_angstrom_squared'] for r in rows),'table_sha256':sha(table),'elapsed_seconds':time.monotonic()-start}
    tmp=rp.with_suffix('.partial');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(rp);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--points',type=int,default=960);a=p.parse_args()
    if a.workers<1 or a.points<100:raise ValueError('Invalid resources/resolution')
    checked_receipt(a.snapshot);models=json.loads((a.snapshot/'model_provenance.json').read_text())
    if len({m['model_id'] for m in models})!=len(models):raise ValueError('Duplicate model identities')
    a.output.mkdir(parents=True,exist_ok=True);lock=(a.output/'.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    source=Path(inspect.getfile(ShrakeRupley))
    c={'snapshot_receipt_sha256':sha(a.snapshot/'receipt.json'),'script_sha256':sha(Path(__file__)),'biopython_version':Bio.__version__,'numpy_version':np.__version__,'SASA_source_sha256':sha(source),'probe_radius_angstrom':1.4,'sphere_points_per_atom':a.points,'atomic_radii_angstrom':dict(ATOMIC_RADII),'workers':a.workers,'models':len(models),'residues':sum(m['length'] for m in models),'interpretation':'Heavy-atom solvent-accessible area for isolated predicted full chains. No experimental exposure, interface, pocket or categorical core assignment. Confidence retained, not used to remove occluding atoms. Inter-domain orientation and prediction uncertainty can affect accessibility.'}
    cp=a.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=c:raise ValueError('Configuration changed')
    cp.write_text(json.dumps(c,indent=2)+'\n');results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(run_one,(m,str(a.output),sha(cp),a.points)) for m in models]
        for f in as_completed(futures):
            r=f.result();results.append(r);print(len(results),'/',len(models),r['model_id'],r['residues'],round(r['elapsed_seconds'],3),flush=True)
    result={'status':'complete_snapshot_predicted_accessibility','config_sha256':sha(cp),'models':len(results),'residues':sum(r['residues'] for r in results),'entry_receipts':{r['model_id']:sha(a.output/(r['model_id']+'.receipt.json')) for r in sorted(results,key=lambda r:r['model_id'])},'interpretation':c['interpretation']}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
