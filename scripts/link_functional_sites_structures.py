#!/usr/bin/env python3
"""Join conserved and nonconserved Pfam site correspondences to one predictor snapshot."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT,sha
from prepare_paired_phylogenetic_inputs import write_table


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def aggregate_sites(rows):
    grouped=defaultdict(list)
    for r in rows:grouped[r['hit_id'],r['profile_match_position']].append(r)
    sites=[]
    identity=['sequence_id','pfam_accession','pfam_name','pfam_type','protein_residue_1based','observed_amino_acid','overlaps_another_ga_hit','hmm_coverage']
    for (hit,pos),group in sorted(grouped.items()):
        first=group[0]
        if any(any(r[k]!=first[k] for k in identity) for r in group):raise ValueError('Conflicting reference projections for one profile site')
        sites.append({k:first[k] for k in identity}|{'site_id':hit+':'+pos,'hit_id':hit,'profile_match_position':pos,
            'expected_amino_acids':';'.join(sorted({r['expected_amino_acid'] for r in group})),
            'reference_accessions':';'.join(sorted({r['reference_uniprot_accession'] for r in group})),
            'pattern_statuses':';'.join(sorted({r['pattern_status'] for r in group})),
            'conserved_candidate':any(r['candidate_active_site']=='True' for r in group)})
    return sites


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['functional-sites','snapshot','encodings','paired','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--source-label',required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable new site join')
    receipts={k:checked_receipt(getattr(a,k)) for k in ['functional_sites','snapshot','encodings','paired']}
    if receipts['functional_sites']['status']!='complete_marker_pfam_active_site_projection' or receipts['encodings']['status']!='complete_native_3di_feature_audit':raise ValueError('Completed site projection and qualified encodings required')
    if receipts['encodings']['mapping_receipt_sha256']!=sha(a.snapshot/'receipt.json') or receipts['paired']['source_receipts']['snapshot']['sha256']!=sha(a.snapshot/'receipt.json') or receipts['paired']['source_receipts']['encodings']['sha256']!=sha(a.encodings/'receipt.json'):raise ValueError('Prediction mapping/encoding/paired provenance differs')
    inputs=ROOT/'data/domains/marker-inputs-v1';checked_receipt(inputs)
    if sha(inputs/'receipt.json')!=receipts['functional_sites']['input_receipt_sha256']:raise ValueError('Different functional protein inputs')
    sites=aggregate_sites(read(a.functional_sites/'site_projection.tsv'));sites_by_seq=defaultdict(list)
    for site in sites:sites_by_seq[site['sequence_id']].append(site)
    protein_links=[r for r in read(inputs/'protein_links.tsv') if r['sequence_id'] in sites_by_seq]
    link_by_key={(r['marker'],r['taxon_id']):r for r in protein_links}
    if len(link_by_key)!=len(protein_links) or set(sites_by_seq)!={r['sequence_id'] for r in protein_links}:raise ValueError('Missing/repeated functional marker links')
    structures={(r['marker'],r['taxon_id']):r for r in read(a.snapshot/'marker_structure_links.tsv')}
    enc={(r['model_id'],r['version']):r for r in read(a.encodings/'model_summary.tsv')}
    wanted={key:{int(s['protein_residue_1based']) for s in sites_by_seq[r['sequence_id']] if s['protein_residue_1based']} for key,r in link_by_key.items()}
    mapped={}
    with gzip.open(a.snapshot/'matrix_to_structure_residues.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=r['marker'],r['taxon_id'];pos=int(r['protein_residue_1based'])
            if key in wanted and pos in wanted[key]:
                k=(*key,pos)
                if k in mapped:raise ValueError('Repeated functional residue mapping')
                model=structures[key]
                if r['model_id']!=model['model_id'] or r['model_version']!=model['model_version']:raise ValueError('Residue/model identity differs')
                mapped[k]=int(r['matrix_column_1based'])
    matrix_path=Path(receipts['paired']['source_receipts']['matrix']['path']);mr=checked_receipt(matrix_path)
    if sha(matrix_path/'receipt.json')!=receipts['snapshot']['matrix_receipt_sha256']:raise ValueError('Different profile matrix')
    matrix={r.id:str(r.seq) for r in SeqIO.parse(matrix_path/'matrix.faa','fasta')}
    paired={}
    for r in read(a.paired/'marker_summary.tsv'):
        if r['status']!='ready_for_inference':continue
        folder=a.paired/r['marker'];paired[r['marker']]={'columns':{int(x['matrix_column_1based']):int(x['paired_column_1based']) for x in read(folder/'columns.tsv')},
            'aa':{x.id:str(x.seq) for x in SeqIO.parse(folder/'aa.faa','fasta')},'states':{x.id:str(x.seq) for x in SeqIO.parse(folder/'3di.faa','fasta')}}
    cache={};rows=[]
    for key,link in sorted(link_by_key.items()):
        model=structures.get(key);data=None
        if model:
            if model['sequence_sha256']!=link['sequence_sha256'] or model['protein_id']!=link['protein_id']:raise ValueError('Functional/model protein differs')
            mk=model['model_id'],model['model_version']
            if mk not in cache:
                er=enc[mk];path=ROOT/er['encoding_path']
                if sha(path)!=er['encoding_sha256']:raise ValueError('Changed functional structure encoding')
                with np.load(path,allow_pickle=False) as d:cache[mk]={k:d[k].copy() for k in d.files}
                cache[mk]['sequence']=str(cache[mk]['sequence']);cache[mk]['states']=str(cache[mk]['states'])
                if hashlib.sha256(cache[mk]['sequence'].encode()).hexdigest()!=link['sequence_sha256']:raise ValueError('Encoding sequence differs')
            data=cache[mk]
        for site in sites_by_seq[link['sequence_id']]:
            row={**site,'source_label':a.source_label,'marker':key[0],'taxon_id':key[1],'protein_id':link['protein_id'],'model_id':model['model_id'] if model else '',
                'model_version':model['model_version'] if model else '', 'model_available':bool(model),'site_status':'no_model_in_snapshot','matrix_column_1based':'','paired_column_1based':'','focal_plddt':'','native_valid':'','native_state':'','feature_min_plddt':'','feature_max_pae':'','joint_feature_confident':False,'observed_in_paired_alignment':False}
            if not site['protein_residue_1based']:row['site_status']='gap_in_functional_profile_alignment'
            elif data is not None:
                pos=int(site['protein_residue_1based']);i=pos-1
                if not 0<=i<len(data['sequence']) or data['sequence'][i]!=site['observed_amino_acid']:raise ValueError('Functional residue differs from structural sequence')
                valid=bool(data['valid'][i]);conf=float(data['ca_plddt'][i]);column=mapped.get((*key,pos))
                row.update(focal_plddt=conf,native_valid=valid,matrix_column_1based=column or '',site_status='mapped_to_marker_matrix' if column else 'model_site_outside_marker_matrix')
                if valid:
                    low=float(data['feature_min_plddt'][i]);pae=float(data['feature_max_pae'][i])
                    if not np.isfinite([low,pae,conf]).all():raise ValueError('Invalid functional feature confidence')
                    row.update(native_state=data['states'][i],feature_min_plddt=low,feature_max_pae=pae,joint_feature_confident=low>=70 and pae<=10)
                if column:
                    if matrix[key[1]][column-1]!=site['observed_amino_acid']:raise ValueError('Functional matrix residue differs')
                    pair=paired.get(key[0]);pc=pair['columns'].get(column) if pair else None
                    if pc and key[1] in pair['aa'] and pair['aa'][key[1]][pc-1]!='?':
                        if not row['joint_feature_confident'] or pair['aa'][key[1]][pc-1]!=site['observed_amino_acid'] or pair['states'][key[1]][pc-1]!=row['native_state']:raise ValueError('Functional paired AA/state correspondence differs')
                        row.update(paired_column_1based=pc,observed_in_paired_alignment=True)
            rows.append(row)
    a.output.mkdir(parents=True);write_table(a.output/'site_structure_links.tsv',rows)
    observed=[r for r in rows if r['observed_in_paired_alignment']]
    result={'status':'complete_functional_site_structure_alignment_join','source_label':a.source_label,'source_receipts':{k:{'path':str(getattr(a,k)),'sha256':sha(getattr(a,k)/'receipt.json')} for k in receipts},'script_sha256':sha(Path(__file__)),
        'unique_profile_sites':len(sites),'taxon_marker_site_rows':len(rows),'site_status_counts':dict(Counter(r['site_status'] for r in rows)),'paired_observed_rows':len(observed),'paired_observed_conserved_candidate_rows':sum(r['conserved_candidate'] for r in observed),'paired_observed_other_correspondence_rows':sum(not r['conserved_candidate'] for r in observed),'paired_observed_taxa':len({r['taxon_id'] for r in observed}),'paired_observed_markers':len({r['marker'] for r in observed}),
        'interpretation':'All projected functional correspondences retained, including nonconserved/gapped patterns and overlap ambiguity. One prediction source per output. Candidate status is separate from structure availability, native feature confidence and actual observation in matched AA/3Di alignments. Repeated homologous sites across taxa are not independent evolutionary events. No catalytic validation, branch effect, selection or causal claim.',
        'artifacts':{'site_structure_links.tsv':sha(a.output/'site_structure_links.tsv')}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
