#!/usr/bin/env python3
"""Compare all annotated Pfam hit spans in audited experimental-control residue masks."""
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.transform import Rotation
from audit_joint_path_uncertainty import checked, rows, sha
from prepare_paired_phylogenetic_inputs import write_table


def align(x, y):
    xc=x-x.mean(0); yc=y-y.mean(0)
    rot,_=Rotation.align_vectors(yc,xc)
    fitted=rot.apply(xc)+y.mean(0)
    rms=float(np.sqrt(np.mean(np.sum((fitted-y)**2,axis=1))))
    independent=SVDSuperimposer(); independent.set(y,x); independent.run()
    if not math.isclose(rms,independent.get_rms(),rel_tol=1e-8,abs_tol=1e-8):
        raise ValueError('Superposition implementations differ')
    return fitted,rms


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['comparisons','audit','crosswalk','mapping','pfam','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipt=checked(a.comparisons); checked(a.crosswalk); checked(a.pfam)
    audit=json.loads(a.audit.read_text())
    if audit['status']!='passed_full_experimental_control_geometry_readback' or audit['comparison_receipt_sha256']!=sha(a.comparisons/'receipt.json'):
        raise ValueError('Complete matching geometry audit required')
    for name in ['crosswalk','mapping']:
        if receipt['source_receipts'][name]!=sha(getattr(a,name)/'receipt.json'):raise ValueError('Lineage mismatch')
    links={(r['entity_id'],r['alphafold_model_id'],r['esmfold_model_id']):r for r in rows(a.crosswalk/'crosswalk.tsv')}
    lengths={r['sequence_sha256']:int(r['sequence_length']) for r in links.values()}
    hits=defaultdict(list)
    for hit in rows(a.pfam/'raw_annotated_hits.tsv'):
        digest=hit['sequence_id'][1:]
        if digest not in lengths:continue
        if int(hit['protein_length'])!=lengths[digest] or not 1<=int(hit['alignment_start'])<=int(hit['alignment_end'])<=lengths[digest]:
            raise ValueError('Pfam sequence length/span mismatch')
        hits[digest].append(hit)
    mr=json.loads((a.mapping/'receipt.json').read_text()); pred_cache={}; observed_cache={}
    root=Path(__file__).resolve().parents[1]
    def prediction(link,label):
        path=root/link[label+'_path']
        if path not in pred_cache:
            if sha(path)!=link[label+'_sha256']:raise ValueError('Changed prediction')
            d=MMCIF2Dict(str(path)); coords={}
            for atom,i,x,y,z in zip(d['_atom_site.label_atom_id'],d['_atom_site.label_seq_id'],d['_atom_site.Cartn_x'],d['_atom_site.Cartn_y'],d['_atom_site.Cartn_z']):
                if atom=='CA':
                    if int(i) in coords:raise ValueError('Duplicate CA')
                    coords[int(i)]=np.asarray([float(x),float(y),float(z)])
            pred_cache[path]=coords
        return pred_cache[path]
    def observed(row):
        entry=row['entry_id']
        if entry not in observed_cache:
            rp=a.mapping/(entry+'.receipt.json'); r=json.loads(rp.read_text());f=a.mapping/(entry+'.residues.tsv.gz')
            if sha(rp)!=mr['entry_receipt_sha256'][entry] or sha(f)!=r['table_sha256']:raise ValueError('Changed mapping')
            grid=defaultdict(dict)
            with gzip.open(f,'rt') as stream:
                for x in csv.DictReader(stream,delimiter='\t'):
                    if entry+'_'+x['entity_id'] not in {key[0] for key in links} or x['CA_status']!='unambiguous_full_occupancy_CA':continue
                    atom=json.loads(x['atom_records_json']);assert len(atom)==1
                    key=(x['entity_id'],x['model_number'],x['label_asym_id'])
                    grid[key][int(x['label_seq_id'])]=np.array([float(atom[0][k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
            observed_cache[entry]=grid
        return observed_cache[entry][row['entity_id'],row['deposited_model'],row['label_asym_id']]
    output=[];excluded=[];superpositions=0;grids=0
    for source in ['comparisons.tsv','exclusions.tsv']:
        for row in rows(a.comparisons/source):
            grids+=1;digest=row['sequence_sha256'];positions=json.loads(row['residue_positions_json'])
            link=links[(row['entry_id']+'_'+row['entity_id'],row['alphafold_model_id'],row['esmfold_model_id'])]
            if digest!=link['sequence_sha256']:raise ValueError('Sequence mismatch')
            maps=[prediction(link,'alphafold'),prediction(link,'esmfold'),observed(row)]
            xyz=[np.array([m[i] for i in positions]) for m in maps]; global_fits={}
            if len(positions)>=3:
                for name,left,right in [('af_experiment',0,2),('esm_experiment',1,2),('af_esm',0,1)]:
                    global_fits[name]=align(xyz[left],xyz[right])[0];superpositions+=1
            for hit in hits[digest]:
                lo,hi=int(hit['alignment_start']),int(hit['alignment_end']);indices=[i for i,v in enumerate(positions) if lo<=v<=hi]
                out={k:row[k] for k in ['sequence_sha256','alphafold_model_id','esmfold_model_id','entry_id','entity_id','deposited_model','label_asym_id','joint_predicted_plddt_cutoff']}
                out.update({k:hit[k] for k in ['hit_id','pfam_accession','pfam_name','pfam_type','hmm_coverage','alignment_start','alignment_end']})
                out.update(matched_residues=len(indices),span_length=hi-lo+1,whole_protein_eligible=source=='comparisons.tsv')
                if len(indices)<20 or len(indices)*2<hi-lo+1:
                    excluded.append(out|{'reason':'fewer_than_20_CA_or_below_half_hit_span'});continue
                for name,left,right in [('af_experiment',0,2),('esm_experiment',1,2),('af_esm',0,1)]:
                    _,rms=align(xyz[left][indices],xyz[right][indices]);superpositions+=1
                    out[name+'_segment_rmsd']=rms
                    out[name+'_after_whole_mask_fit_rmsd']=float(np.sqrt(np.mean(np.sum((global_fits[name][indices]-xyz[right][indices])**2,axis=1))))
                output.append(out)
    if grids!=receipt['accepted_rows']+receipt['excluded_rows']:raise ValueError('Incomplete source grid')
    a.output.mkdir(parents=True)
    for name,data in [('segments.tsv',output),('exclusions.tsv',excluded)]:
        if data:write_table(a.output/name,data)
    result=dict(status='complete_experimental_control_pfam_segment_geometry',source_grids=grids,
        proteins=len(lengths),proteins_with_hits=len(hits),pfam_hits=sum(map(len,hits.values())),
        accepted_segment_rows=len(output),excluded_segment_rows=len(excluded),independent_superposition_checks=superpositions,
        source_receipts={name:sha(getattr(a,name)/'receipt.json') for name in ['comparisons','crosswalk','mapping','pfam']},
        geometry_audit_sha256=sha(a.audit),script_sha256=sha(Path(__file__)),
        interpretation='All raw Pfam alignment spans retained, including overlapping hits, repeats and families; these are not resolved independent structural domains. Same audited joint residue masks, with segment eligibility >=20 CA and >=half hit span, including rows excluded at whole-protein coverage. Independent segment versus whole-mask rigid fits are descriptive and mechanically favor separate fits; no significance, accuracy or domain-motion causality claim. Independent full grid/mask readback remains pending.',
        artifacts={f.name:sha(f) for f in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
