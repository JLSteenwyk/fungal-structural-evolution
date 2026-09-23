#!/usr/bin/env python3
"""Reconstruct all functional-site output rows from projections and qualified arrays."""
import argparse
from collections import defaultdict,Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from catalog_whole_proteome_structures import sha


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--result',type=Path,required=True);a=ap.parse_args()
    receipt=checked_receipt(a.result)
    if receipt['status']!='complete_functional_site_structure_alignment_join':raise ValueError('Incomplete join')
    paths={k:Path(v['path']) for k,v in receipt['source_receipts'].items()}
    for k,p in paths.items():
        checked_receipt(p)
        if sha(p/'receipt.json')!=receipt['source_receipts'][k]['sha256']:raise ValueError('Source binding differs')
    projection=table(paths['functional_sites']/'site_projection.tsv');groups=defaultdict(list)
    for r in projection:groups[r['hit_id']+':'+r['profile_match_position']].append(r)
    static=['sequence_id','pfam_accession','pfam_name','pfam_type','protein_residue_1based','observed_amino_acid','overlaps_another_ga_hit','hmm_coverage','hit_id','profile_match_position']
    sites={};by_sequence=defaultdict(list)
    for key,rs in groups.items():
        first=rs[0]
        if any(any(x[k]!=first[k] for k in static) for x in rs):raise ValueError('Conflicting projection')
        s={k:first[k] for k in static};s.update(site_id=key,conserved_candidate=str(any(x['candidate_active_site']=='True' for x in rs)))
        for target,source in [('expected_amino_acids','expected_amino_acid'),('reference_accessions','reference_uniprot_accession'),('pattern_statuses','pattern_status')]:s[target]=';'.join(sorted({x[source] for x in rs}))
        sites[key]=s;by_sequence[s['sequence_id']].append(key)
    inputs=Path('data/domains/marker-inputs-v1');checked_receipt(inputs)
    source_receipt=json.loads((paths['functional_sites']/'receipt.json').read_text())
    if sha(inputs/'receipt.json')!=source_receipt['input_receipt_sha256']:raise ValueError('Wrong functional input proteins')
    expected={}
    for link in table(inputs/'protein_links.tsv'):
        for site in by_sequence.get(link['sequence_id'],[]):
            key=link['marker'],link['taxon_id'],site
            if key in expected:raise ValueError('Repeated source identity')
            expected[key]=link
    rows=table(a.result/'site_structure_links.tsv');observed_keys={(r['marker'],r['taxon_id'],r['site_id']) for r in rows}
    if len(observed_keys)!=len(rows) or observed_keys!=set(expected):raise ValueError('Full row universe differs')
    structures={(r['marker'],r['taxon_id']):r for r in table(paths['snapshot']/'marker_structure_links.tsv')}
    wanted={(r['marker'],r['taxon_id'],r['protein_residue_1based']) for r in rows if r['protein_residue_1based']}
    positions={}
    with gzip.open(paths['snapshot']/'matrix_to_structure_residues.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=r['marker'],r['taxon_id'],r['protein_residue_1based']
            if key in wanted:
                if key in positions:raise ValueError('Repeated residue mapping')
                positions[key]=r
    arrays={(r['model_id'],r['version']):r for r in table(paths['encodings']/'model_summary.tsv')};cache={};paired={}
    for r in table(paths['paired']/'marker_summary.tsv'):
        if r['status']=='ready_for_inference':
            d=paths['paired']/r['marker'];paired[r['marker']]=({x.id:str(x.seq) for x in SeqIO.parse(d/'aa.faa','fasta')},{x.id:str(x.seq) for x in SeqIO.parse(d/'3di.faa','fasta')},{x['matrix_column_1based']:int(x['paired_column_1based']) for x in table(d/'columns.tsv')})
    pr=json.loads((paths['paired']/'receipt.json').read_text());matrix_path=Path(pr['source_receipts']['matrix']['path']);checked_receipt(matrix_path)
    if sha(matrix_path/'receipt.json')!=pr['source_receipts']['matrix']['sha256']:raise ValueError('Changed matrix')
    matrix={x.id:str(x.seq) for x in SeqIO.parse(matrix_path/'matrix.faa','fasta')}
    numeric={'focal_plddt','feature_min_plddt','feature_max_pae'};observed=[]
    for row in rows:
        marker,taxon,site=row['marker'],row['taxon_id'],row['site_id'];link=expected[marker,taxon,site];s=sites[site];model=structures.get((marker,taxon));pos=s['protein_residue_1based']
        values={**s,'marker':marker,'taxon_id':taxon,'protein_id':link['protein_id'],'source_label':receipt['source_label'],'model_available':str(model is not None),'model_id':model['model_id'] if model else '', 'model_version':model['model_version'] if model else '', 'site_status':'no_model_in_snapshot','matrix_column_1based':'','paired_column_1based':'','focal_plddt':'','native_valid':'','native_state':'','feature_min_plddt':'','feature_max_pae':'','joint_feature_confident':'False','observed_in_paired_alignment':'False'}
        if model and (model['protein_id']!=link['protein_id'] or model['sequence_sha256']!=link['sequence_sha256']):raise ValueError('Wrong model protein')
        if not pos:values['site_status']='gap_in_functional_profile_alignment'
        elif model:
            mk=model['model_id'],model['model_version'];er=arrays[mk]
            if mk not in cache:
                ep=Path(er['encoding_path'])
                if sha(ep)!=er['encoding_sha256']:raise ValueError('Changed array')
                with np.load(ep,allow_pickle=False) as data:cache[mk]={k:data[k].copy() for k in data.files}
            data=cache[mk];sequence=str(data['sequence']);i=int(pos)-1
            if hashlib.sha256(sequence.encode()).hexdigest()!=link['sequence_sha256'] or sequence[i]!=s['observed_amino_acid']:raise ValueError('Wrong source residue')
            valid=bool(data['valid'][i]);mapping=positions.get((marker,taxon,pos));col=mapping['matrix_column_1based'] if mapping else ''
            values.update(native_valid=str(valid),focal_plddt=float(data['ca_plddt'][i]),matrix_column_1based=col,site_status='mapped_to_marker_matrix' if col else 'model_site_outside_marker_matrix')
            if valid:
                low=float(data['feature_min_plddt'][i]);pae=float(data['feature_max_pae'][i])
                if not np.isfinite([low,pae,values['focal_plddt']]).all():raise ValueError('Nonfinite confidence')
                values.update(native_state=str(data['states'])[i],feature_min_plddt=low,feature_max_pae=pae,joint_feature_confident=str(low>=70 and pae<=10))
            if col:
                if mapping['model_id']!=mk[0] or mapping['model_version']!=mk[1] or matrix[taxon][int(col)-1]!=sequence[i]:raise ValueError('Matrix/model identity differs')
                if marker in paired:
                    aa,ss,cols=paired[marker];pc=cols.get(col)
                    if pc and taxon in aa and aa[taxon][pc-1]!='?':
                        if values['joint_feature_confident']!='True' or aa[taxon][pc-1]!=sequence[i] or ss[taxon][pc-1]!=values['native_state']:raise ValueError('Paired observation differs')
                        values.update(paired_column_1based=str(pc),observed_in_paired_alignment='True')
        if set(values)!=set(row):raise ValueError('Unchecked output fields')
        for field,value in values.items():
            if field in numeric and value!='':
                if not math.isclose(float(row[field]),value,rel_tol=1e-12,abs_tol=1e-12):raise ValueError('Numeric mismatch: '+field)
            elif row[field]!=str(value):raise ValueError('Field mismatch: '+field)
        if values['observed_in_paired_alignment']=='True':observed.append(row)
    counts={'taxon_marker_site_rows':len(rows),'unique_profile_sites':len(sites),'paired_observed_rows':len(observed),'paired_observed_conserved_candidate_rows':sum(r['conserved_candidate']=='True' for r in observed),'paired_observed_other_correspondence_rows':sum(r['conserved_candidate']!='True' for r in observed),'paired_observed_taxa':len({r['taxon_id'] for r in observed}),'paired_observed_markers':len({r['marker'] for r in observed}),'site_status_counts':dict(Counter(r['site_status'] for r in rows))}
    if any(receipt[k]!=v for k,v in counts.items()):raise ValueError('Summary differs')
    result={'status':'passed_full_completed_functional_site_row_readback',**counts,'producer_receipt_sha256':sha(a.result/'receipt.json'),'script_sha256':sha(__file__),'scope':'Every output field and full row universe reconstructed from original projected correspondences, protein links, mapped positions, qualified arrays and paired FASTA. No producer join functions imported. Shared prior projections, array qualification, file parsers and hash helper; no new coordinate, PAE, HMM or biological validation.'}
    (a.result/'readback.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
