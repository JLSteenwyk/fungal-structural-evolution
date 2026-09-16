#!/usr/bin/env python3
"""Independently check every Pfam segment grid and summarize nested geometry medians."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['segments','comparisons','pfam','output','summary']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.summary.exists():raise FileExistsError('Use fresh outputs')
    for folder in [a.segments,a.comparisons,a.pfam]:
        r=json.loads((folder/'receipt.json').read_text())
        for name,h in r['artifacts'].items():
            if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=h:raise ValueError('Changed source artifact')
    read=lambda f:pd.read_csv(f,sep='\t',dtype={'deposited_model':str,'entity_id':str})
    accepted=read(a.segments/'segments.tsv');excluded=read(a.segments/'exclusions.tsv')
    source=pd.concat([read(a.comparisons/f) for f in ['comparisons.tsv','exclusions.tsv']],ignore_index=True)
    hits=read(a.pfam/'raw_annotated_hits.tsv');hits=hits[hits.sequence_id.isin('S'+source.sequence_sha256)].copy()
    hits['sequence_sha256']=hits.sequence_id.str[1:]
    keys=['entry_id','entity_id','deposited_model','label_asym_id','alphafold_model_id','esmfold_model_id','joint_predicted_plddt_cutoff','hit_id'];expected={}
    for r in source.merge(hits,on='sequence_sha256').to_dict('records'):
        n=sum(int(r['alignment_start'])<=i<=int(r['alignment_end']) for i in json.loads(r['residue_positions_json']))
        width=int(r['alignment_end'])-int(r['alignment_start'])+1;key=tuple(r[k] for k in keys)
        if key in expected:raise ValueError('Duplicate source grid')
        expected[key]=(n,width,n>=20 and n*2>=width)
    seen=set()
    for status,df in [(True,accepted),(False,excluded)]:
        for r in df.to_dict('records'):
            key=tuple(r[k] for k in keys)
            if key in seen or expected[key]!=(r['matched_residues'],r['span_length'],status):raise ValueError('Grid count/eligibility differs')
            seen.add(key)
    if seen!=expected.keys():raise ValueError('Incomplete grid')
    metrics=[c for c in accepted if c.endswith('_rmsd')]
    group=['sequence_sha256','alphafold_model_id','hit_id','pfam_name','pfam_type','alignment_start','alignment_end','joint_predicted_plddt_cutoff']
    x=accepted
    for keys in [group+['entry_id','deposited_model'],group+['entry_id'],group]:
        x=x.groupby(keys,dropna=False)[metrics].median().reset_index()
    a.summary.parent.mkdir(parents=True,exist_ok=True);x.to_csv(a.summary,sep='\t',index=False)
    r=dict(status='passed_full_segment_grid_count_eligibility_readback',source_threshold_grids=len(source),
        pfam_hits=len(hits),segment_grids=len(expected),accepted_rows=len(accepted),excluded_rows=len(excluded),
        method='Independent pandas join of all source threshold rows to sequence-keyed Pfam hits; exact identity universe, boundary-derived matched counts, span lengths and eligibility reconstructed for every row. Geometry independently checked by two superposition implementations in producer; no independent full geometry replay here.',
        source_receipt_sha256=hashlib.sha256((a.segments/'receipt.json').read_bytes()).hexdigest(),
        summary_sha256=hashlib.sha256(a.summary.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
