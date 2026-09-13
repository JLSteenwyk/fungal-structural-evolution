#!/usr/bin/env python3
"""Audit paired resolution artifacts and quantify zero-area classification sensitivity."""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
import numpy as np
from audit_busco_gene_copies import sha,read_table
from assess_pae_sensitivity import checked_receipt
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--assessment',type=Path,required=True);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable audit')
    result=checked_receipt(a.assessment);checked_receipt(a.snapshot);config=json.loads((a.assessment/'config.json').read_text())
    if result['config_sha256']!=sha(a.assessment/'config.json') or config['snapshot_receipt_sha256']!=sha(a.snapshot/'receipt.json'):raise ValueError('Source lineage mismatch')
    models=json.loads((a.snapshot/'model_provenance.json').read_text());expected=[]
    for low,high in config['length_bins']:
        candidates=[m for m in models if low<=m['length']<=high]
        expected.extend(sorted(candidates,key=lambda m:hashlib.sha256(('accessibility_resolution_v1|'+m['model_id']).encode()).hexdigest())[:2])
    if expected!=config['selected_models']:raise ValueError('Selection differs')
    stats={r['model_id']:r for r in read_table(a.assessment/'model_resolution_summary.tsv')};rows=[];all_delta=[];pin_count=0
    if set(stats)!={m['model_id'] for m in expected}:raise ValueError('Summary model universe differs')
    for m in expected:
        sid=m['model_id'];vectors=[];confidence=[]
        for points in [960,3840]:
            folder=a.assessment/str(points);cp=folder/'config.json';c=json.loads(cp.read_text());rp=folder/(sid+'.receipt.json');r=json.loads(rp.read_text());table=folder/(sid+'.residues.tsv.gz')
            if c!={'parent_config_sha256':sha(a.assessment/'config.json'),'points':points}:raise ValueError('Resolution config differs')
            if result['comparison_receipt_sha256'][str(folder/sid)]!=sha(rp) or r['config_sha256']!=sha(cp) or r['table_sha256']!=sha(table) or r['model_sha256']!=m['sha256']:raise ValueError('Artifact lineage differs')
            pin_count+=1
            with gzip.open(table,'rt') as f:values=list(csv.DictReader(f,delimiter='\t'))
            if len(values)!=m['length'] or [int(x['protein_residue_1based']) for x in values]!=list(range(1,m['length']+1)):raise ValueError('Residue grid differs')
            if hashlib.sha256(''.join(x['amino_acid'] for x in values).encode()).hexdigest()!=m['sequence_sha256']:raise ValueError('Sequence differs')
            v=np.array([float(x['sasa_angstrom_squared']) for x in values]);plddt=np.array([float(x['ca_plddt']) for x in values])
            if not np.isfinite(v).all() or (v<0).any() or not math.isclose(v.sum(),r['total_sasa_angstrom_squared'],rel_tol=1e-12):raise ValueError('Invalid ASA values/totals')
            vectors.append(v);confidence.append(plddt)
        if not np.array_equal(*confidence):raise ValueError('Prediction confidence changed across resolution')
        x,y=vectors;delta=abs(x-y);all_delta.extend(delta.tolist());s=stats[sid]
        computed={'median_absolute_residue_ASA_difference_A2':np.median(delta),'p95_absolute_residue_ASA_difference_A2':np.quantile(delta,.95),'maximum_absolute_residue_ASA_difference_A2':delta.max(),'total_ASA_960_A2':x.sum(),'total_ASA_3840_A2':y.sum()}
        if any(not math.isclose(v,float(s[k]),rel_tol=1e-12,abs_tol=1e-10) for k,v in computed.items()):raise ValueError('Reported summary differs')
        changed=(x==0)!=(y==0);high=confidence[0]>=70
        rows.append({'model_id':sid,'length':len(x),'zero_ASA_960':int((x==0).sum()),'zero_ASA_3840':int((y==0).sum()),'zero_nonzero_disagreements':int(changed.sum()),'high_confidence_residues':int(high.sum()),'high_confidence_zero_nonzero_disagreements':int((changed&high).sum()),'maximum_difference_for_zero_nonzero_disagreements_A2':float(delta[changed].max()) if changed.any() else 0,**{k:float(v) for k,v in computed.items()}})
    if pin_count!=len(result['comparison_receipt_sha256']) or len(rows)!=result['models'] or len(all_delta)!=result['residues']:raise ValueError('Completion counts differ')
    a.output.mkdir(parents=True);write_table(a.output/'resolution_sensitivity.tsv',rows)
    r={'status':'passed_all_selected_accessibility_resolution_outputs','assessment_receipt_sha256':sha(a.assessment/'receipt.json'),'snapshot_receipt_sha256':sha(a.snapshot/'receipt.json'),'script_sha256':sha(Path(__file__)),'models':len(rows),'residues':len(all_delta),'median_absolute_residue_difference_A2':float(np.median(all_delta)),'p95_absolute_residue_difference_A2':float(np.quantile(all_delta,.95)),'maximum_absolute_residue_difference_A2':max(all_delta),'zero_nonzero_disagreements':sum(r['zero_nonzero_disagreements'] for r in rows),'interpretation':'All selected input identities, output hashes/grids/sequences and reported resolution summaries independently checked. Zero/nonzero area sensitivity is a numerical diagnostic, not core/surface biological truth. Higher resolution is not exact ground truth; selected models do not prove full-dataset convergence.','artifacts':{'resolution_sensitivity.tsv':sha(a.output/'resolution_sensitivity.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
