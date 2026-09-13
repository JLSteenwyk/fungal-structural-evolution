#!/usr/bin/env python3
"""Assess surface sampling resolution with fixed length-stratified model selection."""
import argparse,csv,gzip,hashlib,json,inspect
from pathlib import Path
from collections import defaultdict
import Bio,numpy as np
from Bio.PDB.SASA import ShrakeRupley,ATOMIC_RADII
from annotate_predicted_accessibility import run_one
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use new immutable resolution assessment')
    checked_receipt(a.snapshot);models=json.loads((a.snapshot/'model_provenance.json').read_text());bins=[(1,128),(129,256),(257,512),(513,1024),(1025,1000000)];selected=[]
    for low,high in bins:
        candidates=[m for m in models if low<=m['length']<=high]
        selected.extend(sorted(candidates,key=lambda m:hashlib.sha256(('accessibility_resolution_v1|'+m['model_id']).encode()).hexdigest())[:2])
    a.output.mkdir(parents=True);config={'snapshot_receipt_sha256':sha(a.snapshot/'receipt.json'),'script_sha256':sha(Path(__file__)),'producer_sha256':sha(Path(inspect.getfile(run_one))),'SASA_source_sha256':sha(Path(inspect.getfile(ShrakeRupley))),'biopython_version':Bio.__version__,'numpy_version':np.__version__,'atomic_radii_angstrom':dict(ATOMIC_RADII),'probe_radius_angstrom':1.4,'sphere_points':[960,3840],'selection':'Two smallest SHA256(accessibility_resolution_v1|model_id) models per fixed length bin, independent of ASA values; all bins retained when nonempty.','length_bins':bins,'selected_models':selected}
    cp=a.output/'config.json';cp.write_text(json.dumps(config,indent=2)+'\n');summary=[];pins={}
    for model in selected:
        values=[]
        for points in [960,3840]:
            folder=a.output/str(points);folder.mkdir(exist_ok=True)
            # Each resolution has a distinct configuration hash, including point count.
            cfg=folder/'config.json';cfg.write_text(json.dumps({'parent_config_sha256':sha(cp),'points':points},indent=2)+'\n')
            r=run_one((model,str(folder),sha(cfg),points));pins[str(folder/model['model_id'])]=sha(folder/(model['model_id']+'.receipt.json'))
            with gzip.open(folder/(model['model_id']+'.residues.tsv.gz'),'rt') as f:values.append(list(csv.DictReader(f,delimiter='\t')))
        low,high=values
        if [(r['protein_residue_1based'],r['amino_acid']) for r in low]!=[(r['protein_residue_1based'],r['amino_acid']) for r in high]:raise ValueError('Resolution comparison identity mismatch')
        delta=np.abs(np.array([float(r['sasa_angstrom_squared']) for r in low])-np.array([float(r['sasa_angstrom_squared']) for r in high]))
        summary.append({'model_id':model['model_id'],'length':model['length'],'median_absolute_residue_ASA_difference_A2':float(np.median(delta)),'p95_absolute_residue_ASA_difference_A2':float(np.quantile(delta,.95)),'maximum_absolute_residue_ASA_difference_A2':float(delta.max()),'total_ASA_960_A2':sum(float(r['sasa_angstrom_squared']) for r in low),'total_ASA_3840_A2':sum(float(r['sasa_angstrom_squared']) for r in high)})
        print(model['model_id'],summary[-1],flush=True)
    write_table(a.output/'model_resolution_summary.tsv',summary)
    result={'status':'complete_sampled_accessibility_resolution_assessment','config_sha256':sha(cp),'models':len(summary),'residues':sum(r['length'] for r in summary),'comparison_receipt_sha256':pins,'interpretation':'Numerical sampling sensitivity on predetermined length-stratified examples; 3840 points are a higher-resolution comparator, not exact ground truth. Does not establish convergence for every production model, remove prediction error or validate biological exposure.','artifacts':{'model_resolution_summary.tsv':sha(a.output/'model_resolution_summary.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
