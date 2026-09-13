#!/usr/bin/env python3
"""Independently reconstruct every local marker residue mapping from source alignments."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from Bio import AlignIO, SeqIO
from audit_busco_gene_copies import ROOT, sha, read_table


def checked(folder):
    r=json.loads((folder/'receipt.json').read_text())
    for name,digest in r['artifacts'].items():
        if sha(folder/name)!=digest:raise ValueError('Changed artifact: '+name)
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['mapping','prediction-audit','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--expected-links',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable residue audit')
    mapping=checked(a.mapping);audit=checked(a.prediction_audit)
    if mapping['source_policy']['provider']!='local' or mapping['source_policy']['tool']!='ESMFold v1':
        raise ValueError('Local ESMFold source required')
    matrix=ROOT/'results/phylogeny/profile-matrix-50-v1';mr=checked(matrix)
    if sha(matrix/'receipt.json')!=mapping['matrix_receipt_sha256']:raise ValueError('Matrix receipt differs')
    with (matrix/'matrix.faa').open() as handle:
        matrix_sequences={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    models=json.loads((a.mapping/'model_provenance.json').read_text())
    audited={r['sequence_id']:r for r in read_table(a.prediction_audit/'predictions.tsv')}
    model_data={};sequences=set()
    for model in models:
        key=(model['model_id'],str(model['version']))
        if key in model_data:raise ValueError('Duplicate model identity')
        sid='S'+model['sequence_sha256'];sequences.add(sid)
        if sid not in audited or model['prediction_receipt_sha256']!=audited[sid]['prediction_receipt_sha256']:
            raise ValueError('Model outside independent prediction audit')
        for path_key,hash_key in [('path','sha256'),('prediction_receipt_path','prediction_receipt_sha256'),('local_pae_npz_path','local_pae_npz_sha256')]:
            if sha(ROOT/model[path_key])!=model[hash_key]:raise ValueError('Changed model or prediction source')
        with np.load(ROOT/model['local_pae_npz_path'],allow_pickle=False) as data:
            sequence=str(data['sequence']);confidence=data['ca_plddt'].copy()
        if hashlib.sha256(sequence.encode()).hexdigest()!=model['sequence_sha256'] or len(sequence)!=model['length']:
            raise ValueError('NPZ sequence identity differs')
        if confidence.shape!=(len(sequence),) or not np.isfinite(confidence).all() or np.any((confidence<0)|(confidence>100)):
            raise ValueError('Invalid NPZ confidence')
        model_data[key]=(sequence,confidence)
    if sequences!=set(audited) or len(models)!=mapping['distinct_models'] or len(models)!=audit['predictions']:
        raise ValueError('Full local model universe differs')
    links=read_table(a.mapping/'marker_structure_links.tsv');original=read_table(a.prediction_audit/'taxon_links.tsv')
    fields=['marker','taxon_id','protein_id','sequence_sha256'];key=lambda r:tuple(r[k] for k in fields)
    expected_links=read_table(a.expected_links)
    source_links_path=ROOT/'results/phylogeny/markers-full-v1/protein_mapping.tsv'
    source_links=[r for r in read_table(source_links_path) if 'S'+r['sequence_sha256'] in sequences]
    observed=Counter(map(key,links));declared=Counter(map(key,expected_links))
    if (observed!=declared or declared!=Counter(map(key,source_links))
            or Counter(map(key,original))-observed):
        raise ValueError('Full exact-sequence link grid or originating subset differs')
    grouped=defaultdict(list);lookup={}
    for row in links:
        k=(row['marker'],row['taxon_id'])
        if k in lookup:raise ValueError('Repeated marker/taxon link')
        lookup[k]=row;grouped[row['marker']].append(row)
    sites=defaultdict(list)
    for row in read_table(matrix/'site_mapping.tsv'):
        sites[row['marker']].append((int(row['alignment_column_1based']),int(row['matrix_column_1based'])))
    collection=ROOT/'results/phylogeny/profile-alignments-full-v1/receipt.json'
    profile_rows={r['marker']:r for r in json.loads(collection.read_text())['alignments']}
    expected={};input_hashes={str(collection.relative_to(ROOT)):sha(collection)}
    for marker,marker_links in sorted(grouped.items()):
        prefix=ROOT/'results/phylogeny/profile-alignments-full-v1'/marker
        pr=json.loads(prefix.with_suffix('.receipt.json').read_text())
        if pr!=profile_rows[marker] or sha(prefix.with_suffix('.sto'))!=pr['stockholm_sha256']:
            raise ValueError('Profile collection or Stockholm input changed')
        input_hashes[str(prefix.with_suffix('.receipt.json').relative_to(ROOT))]=sha(prefix.with_suffix('.receipt.json'))
        with prefix.with_suffix('.sto').open() as handle:
            aligned={r.id:str(r.seq) for r in AlignIO.read(handle,'stockholm')}
        for link in marker_links:
            sequence,confidence=model_data[(link['model_id'],link['model_version'])]
            text=aligned[link['taxon_id']]
            present=np.array([c not in '-.' for c in text],dtype=bool)
            if ''.join(c for c in text if c not in '-.').upper()!=sequence:
                raise ValueError('Aligned source sequence differs from NPZ')
            positions=np.cumsum(present)
            for profile_col,matrix_col in sites[marker]:
                stock_col=pr['retained_stockholm_columns_1based'][profile_col-1]-1
                if not present[stock_col]:continue
                position=int(positions[stock_col]);k=(marker,link['taxon_id'],matrix_col)
                if k in expected:raise ValueError('Duplicate expected residue identity')
                if matrix_sequences[link['taxon_id']][matrix_col-1]!=sequence[position-1]:
                    raise ValueError('Matrix amino acid differs from mapped source residue')
                expected[k]=(position,float(confidence[position-1]))
    expected_count=len(expected);summaries=defaultdict(list);seen=0
    with gzip.open(a.mapping/'matrix_to_structure_residues.tsv.gz','rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            k=(row['marker'],row['taxon_id'],int(row['matrix_column_1based']))
            if k not in expected:raise ValueError('Unexpected or repeated residue identity')
            position,confidence=expected.pop(k);link=lookup[k[:2]]
            value=float(row['ca_plddt'])
            if (int(row['protein_residue_1based'])!=position or not math.isfinite(value)
                    or abs(value-confidence)>.0051 or row['protein_id']!=link['protein_id']
                    or row['model_id']!=link['model_id'] or row['model_version']!=link['model_version']):
                raise ValueError('Mapped residue position/confidence/source differs')
            summaries[k[:2]].append(value);seen+=1
    if expected or seen!=mapping['matrix_residue_links']:raise ValueError('Missing mapped residues or count mismatch')
    for k,link in lookup.items():
        values=summaries[k]
        if len(values)!=int(link['retained_marker_residues']):raise ValueError('Per-link residue count differs')
        if values:
            if abs(math.fsum(values)/len(values)-float(link['mean_retained_ca_plddt']))>1e-8:
                raise ValueError('Per-link confidence mean differs')
            if abs(sum(v>=70 for v in values)/len(values)-float(link['fraction_retained_ca_plddt_ge70']))>1e-12:
                raise ValueError('Per-link confidence fraction differs')
    a.output.mkdir(parents=True)
    (a.output/'profile_input_hashes.json').write_text(json.dumps(input_hashes,indent=2)+'\n')
    r={'status':'passed_full_local_residue_mapping_readback','models':len(models),'marker_links':len(links),
       'matrix_residue_links':seen,'reconstructed_expected_residue_links':expected_count,
       'mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'prediction_audit_receipt_sha256':sha(a.prediction_audit/'receipt.json'),
       'expected_links_sha256':sha(a.expected_links),'global_marker_links_sha256':sha(source_links_path),
       'originating_links':len(original),'additional_exact_sequence_links':len(links)-len(original),
       'script_sha256':sha(Path(__file__)),
       'scope':'All source/link identities and NPZ/coordinate hashes checked. Every expected non-gap retained Stockholm position independently reconstructed using cumulative residue counts; all exported positions, matrix amino acids and NPZ confidence checked with PDB rounding tolerance. Shared Biopython alignment parser; no independent alignment inference. PAE/native context and biological interpretation remain separate.',
       'artifacts':{p.name:sha(p) for p in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))


if __name__=='__main__':main()
