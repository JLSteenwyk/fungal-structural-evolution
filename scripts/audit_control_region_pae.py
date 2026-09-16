#!/usr/bin/env python3
"""Independently reconstruct every region-PAE mask and quantile using scalar lists."""
import argparse
import csv
import json
import statistics
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from audit_joint_path_uncertainty import checked, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['summary','controls','crosswalk','pfam','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    receipt=checked(a.summary)
    for n in ['controls','crosswalk','pfam']:
        checked(getattr(a,n))
        if receipt['source_receipts'][n]!=sha(getattr(a,n)/'receipt.json'):raise ValueError('Source lineage differs')
    def rows(f):return list(csv.DictReader(f.open(),delimiter='\t'))
    data=rows(a.summary/'region_pae.tsv');models={x['sequence_sha256']:x for x in json.loads((a.controls/'model_provenance.json').read_text())}
    links={x['sequence_sha256']:x for x in rows(a.crosswalk/'crosswalk.tsv')};hits={x['hit_id']:x for x in rows(a.pfam/'raw_annotated_hits.tsv')}
    cache={};quantiles=0;pairvalues=0;parser=MMCIFParser(QUIET=True)
    for row in data:
        sid=row['sequence_sha256']
        if sid not in cache:
            model=models[sid];file=Path(model['local_pae_npz_path']);afpath=Path(links[sid]['alphafold_path'])
            if sha(file)!=model['local_pae_npz_sha256'] or sha(afpath)!=links[sid]['alphafold_sha256']:raise ValueError('Changed model artifact')
            with np.load(file,allow_pickle=False) as d:pae=d['pae'].copy();esm=d['ca_plddt'].copy()
            af={x.id[1]:x['CA'].bfactor for x in parser.get_structure('af',afpath).get_residues()};cache[sid]=(pae,esm,af)
        pae,esm,af=cache[sid];left=hits[row['left_hit_id']];right=hits[row['right_hit_id']];cut=int(row['joint_predicted_plddt_cutoff'])
        spans=[list(range(int(h['alignment_start']),int(h['alignment_end'])+1)) for h in [left,right]]
        selected=[[i-1 for i in s if float(esm[i-1])>=cut and af[i]>=cut] for s in spans]
        within=left['hit_id']==right['hit_id'];overlap=bool(set(spans[0])&set(spans[1])) and not within
        status='excluded_overlapping_hits' if overlap else 'excluded_coverage' if any(len(z)<20 or len(z)*2<len(s) for z,s in zip(selected,spans)) else 'included'
        if status!=row['status'] or [len(s) for s in selected]!=[int(row['left_residues']),int(row['right_residues'])]:raise ValueError('Mask/disposition differs')
        for name,first,second in [('row_left_column_right',selected[0],selected[1]),('row_right_column_left',selected[1],selected[0])]:
            if status!='included':
                if int(row[name+'_pairs'])!=0:raise ValueError('Excluded pair count')
                continue
            values=[float(pae[i,j]) for i in first for j in second if not within or i!=j]
            if len(values)!=int(row[name+'_pairs']):raise ValueError('Pair count differs')
            pairvalues+=len(values)
            np.testing.assert_allclose([float(row[name+'_median_angstrom']),float(row[name+'_q90_angstrom'])],
                [statistics.median(values),statistics.quantiles(values,n=10,method='inclusive')[8]],rtol=1e-12,atol=1e-12)
            quantiles+=2
    r=dict(status='passed_all_control_pae_summary_rows',rows=len(data),quantiles_checked=quantiles,directed_pair_values_checked=pairvalues,
        source_receipt_sha256=sha(a.summary/'receipt.json'),script_sha256=sha(Path(__file__)),
        method='Independent Bio.PDB AF-confidence parsing, Python scalar residue masking and explicit matrix-value lists; statistics median/inclusive quantiles checked against float64 NumPy summaries. Every emitted row disposition and pair count checked. Does not independently establish completeness of emitted region-pair universe or PAE calibration.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
