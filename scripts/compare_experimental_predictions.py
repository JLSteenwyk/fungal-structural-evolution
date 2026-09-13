#!/usr/bin/env python3
"""Descriptive CA agreement with exact-sequence experimental chains, preserving every model."""
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import pdist
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from compare_marker_structures import geometry
from prepare_paired_phylogenetic_inputs import write_table
THREE={k.upper():v for k,v in protein_letters_3to1.items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['mapping','screen','predictions','references','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable comparison output')
    for name in ['screen','predictions','references']:checked_receipt(getattr(a,name))
    mr=json.loads((a.mapping/'receipt.json').read_text());mc=json.loads((a.mapping/'config.json').read_text())
    if mr['status']!='complete_experimental_CA_correspondence_screen' or mr['config_sha256']!=sha(a.mapping/'config.json') or mc['screen_receipt_sha256']!=sha(a.screen/'receipt.json'):raise ValueError('Mapping lineage mismatch')
    sr=json.loads((a.screen/'receipt.json').read_text())
    if sr['mapping_receipt_sha256']!=sha(a.predictions/'receipt.json'):raise ValueError('Prediction snapshot mismatch')
    models={r['model_id']:r for r in json.loads((a.predictions/'model_provenance.json').read_text())}
    refs={r['entry_id']:r for r in read_table(a.references/'entries.tsv')}
    ref_entities={(r['entity_id'],r['model_id']):r for r in read_table(a.references/'entities.tsv')}
    targets=defaultdict(list)
    for r in read_table(a.screen/'sequence_correspondence.tsv'):
        if r['sequence_class']=='exact_full_sequence':
            if ref_entities[r['entity_id'],r['model_id']]['model_sequence_sha256']!=r['model_sequence_sha256']:raise ValueError('Reference target differs')
            targets[r['entity_id']].append(r)
    cache={};results=[];exclusions=[]
    def predicted(target):
        name=target['model_id'];model=models[name]
        if model['sequence_sha256']!=target['model_sequence_sha256']:raise ValueError('Prediction sequence hash mismatch')
        if name not in cache:
            path=ROOT/model['path']
            if sha(path)!=model['sha256']:raise ValueError('Changed prediction file')
            d=MMCIF2Dict(str(path));atoms={}
            for atom,pos,aa,x,y,z,b in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.label_comp_id'],d['_atom_site.Cartn_x'],d['_atom_site.Cartn_y'],d['_atom_site.Cartn_z'],d['_atom_site.B_iso_or_equiv']):
                if atom=='CA':
                    pos=int(pos)
                    if pos in atoms:raise ValueError('Ambiguous predicted CA')
                    atoms[pos]=(THREE.get(aa,'?'),np.array([float(x),float(y),float(z)]),float(b))
            if sorted(atoms)!=list(range(1,model['length']+1)):raise ValueError('Incomplete predicted CA sequence')
            if hashlib.sha256(''.join(atoms[i][0] for i in sorted(atoms)).encode()).hexdigest()!=model['sequence_sha256']:raise ValueError('Predicted CA sequence differs')
            if any(not np.isfinite(v[1]).all() or not math.isfinite(v[2]) or not 0<=v[2]<=100 for v in atoms.values()):raise ValueError('Invalid predicted coordinate/confidence')
            cache[name]=atoms
        return cache[name]
    for entry,pin in sorted(mr['entry_receipt_sha256'].items()):
        rp=a.mapping/(entry+'.receipt.json');r=json.loads(rp.read_text());table=a.mapping/(entry+'.residues.tsv.gz')
        if sha(rp)!=pin or r['config_sha256']!=mr['config_sha256'] or sha(table)!=r['table_sha256']:raise ValueError('Changed mapping table')
        if refs[entry]['methodology']!='experimental':raise ValueError('Nonexperimental reference requires separate analysis')
        groups=defaultdict(list)
        with gzip.open(table,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):groups[row['entity_id'],row['model_number'],row['label_asym_id']].append(row)
        for (entity,number,chain),rows in sorted(groups.items()):
            rows.sort(key=lambda r:int(r['label_seq_id']));positions=[int(r['label_seq_id']) for r in rows]
            if positions!=list(range(1,len(rows)+1)):raise ValueError('Incomplete/duplicate experimental grid')
            sequence=''.join(r['target_aa'] for r in rows)
            for target in targets[entry+'_'+entity]:
                if hashlib.sha256(sequence.encode()).hexdigest()!=target['model_sequence_sha256']:raise ValueError('Experimental grid sequence differs')
                pred=predicted(target);observed={}
                for row in rows:
                    if row['CA_status']=='unambiguous_full_occupancy_CA':
                        atoms=json.loads(row['atom_records_json'])
                        if len(atoms)!=1:raise ValueError('Ambiguous accepted CA')
                        atom=atoms[0];observed[int(row['label_seq_id'])]=np.array([float(atom[k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
                for cutoff in [0,70,90]:
                    selected=sorted(i for i in observed if pred[i][2]>=cutoff);n=len(selected)
                    row={'entry_id':entry,'entity_id':entity,'deposited_model':number,'label_asym_id':chain,'predicted_model_id':target['model_id'],'sequence_sha256':target['model_sequence_sha256'],'sequence_length':len(rows),'experimental_unambiguous_CA':len(observed),'predicted_plddt_cutoff':cutoff,'matched_residues':n,'fraction_full_sequence':n/len(rows),'experimental_method':refs[entry]['methods'],'resolution_combined_json':refs[entry]['resolution_combined_json'],'initial_release_date':refs[entry]['initial_release_date'],'polymer_composition':refs[entry]['polymer_composition'],'interpretation':'descriptive_agreement_pending_experimental_quality_context_and_training_review'}
                    if n<50 or 2*n<len(rows):exclusions.append(row|{'reason':'fewer_than_50_CA_or_below_half_full_sequence'});continue
                    x=np.array([pred[i][1] for i in selected]);y=np.array([observed[i] for i in selected]);pos=np.array(selected);row.update(geometry(x,y,pos,pos))
                    xc=x-x.mean(0);yc=y-y.mean(0);rot,_=Rotation.align_vectors(yc,xc);rms=np.sqrt(np.mean(np.sum((rot.apply(xc)-yc)**2,axis=1)))
                    dx,dy=pdist(x),pdist(y);ii,jj=np.triu_indices(n,1);local=((dx<=15)|(dy<=15))&(abs(pos[ii]-pos[jj])>=3);delta=abs(dx[local]-dy[local])
                    check={'ca_superposition_rmsd_angstrom':rms,'local_distance_pairs':len(delta)}
                    if len(delta):check.update(local_distance_mean_absolute_change_angstrom=delta.mean(),local_distance_rms_change_angstrom=np.sqrt(np.mean(delta**2)))
                    if any(not math.isclose(float(row[k]),v,rel_tol=1e-8,abs_tol=1e-8) for k,v in check.items()):raise ValueError('Independent geometry mismatch')
                    results.append(row)
        print(entry,'compared',flush=True)
    a.output.mkdir(parents=True);write_table(a.output/'comparisons.tsv',results)
    if exclusions:write_table(a.output/'exclusions.tsv',exclusions)
    summary=[]
    for cutoff in [0,70,90]:
        rs=[r for r in results if r['predicted_plddt_cutoff']==cutoff]
        summary.append({'predicted_plddt_cutoff':cutoff,'accepted_chain_model_comparisons':len(rs),'excluded_chain_model_comparisons':sum(r['predicted_plddt_cutoff']==cutoff for r in exclusions),'distinct_predicted_proteins':len({r['predicted_model_id'] for r in rs}),'entries':len({r['entry_id'] for r in rs}),'comparison_weighted_median_CA_RMSD_angstrom':float(np.median([r['ca_superposition_rmsd_angstrom'] for r in rs])) if rs else None})
    receipt={'status':'complete_descriptive_experimental_prediction_agreement','source_receipts':{n:sha(getattr(a,n)/'receipt.json') for n in ['mapping','screen','predictions','references']},'script_sha256':sha(Path(__file__)),'geometry_helper_sha256':sha(Path(__file__).with_name('compare_marker_structures.py')),'accepted_rows':len(results),'excluded_rows':len(exclusions),'threshold_summary':summary,'independent_geometry_recalculations':len(results),'interpretation':'Every experimental chain and deposited model retained; no best-agreement selection. Experimental B factors never used as pLDDT. Only predicted focal pLDDT thresholds, not six-residue or PAE filters. Coverage requires 50 CA and half full sequence. Chain/model-weighted summaries are descriptive and nonindependent; changing thresholds changes residues and sometimes cohort. Partial reference sampling, experimental quality/context, refinement/template overlap and prediction training independence remain unresolved. Agreement is not unbiased prediction accuracy or an evolutionary change estimate.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
