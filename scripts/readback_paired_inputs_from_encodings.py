#!/usr/bin/env python3
"""Reconstruct every eligible paired character independently from qualified arrays."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import ROOT, sha

ALPHABET=set('ACDEFGHIKLMNPQRSTVWY')


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))


def fasta(path):
    rows=list(SeqIO.parse(path,'fasta'))
    result={r.id:str(r.seq) for r in rows}
    if len(result)!=len(rows):
        raise ValueError('Repeated FASTA taxon')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    produced=checked_receipt(a.inputs)
    if produced['status']!='complete_paired_phylogenetic_input_preparation':
        raise ValueError('Paired preparation incomplete')
    sources={k:Path(v['path']) for k,v in produced['source_receipts'].items()}
    pins={str(a.inputs/'receipt.json'):sha(a.inputs/'receipt.json')}
    for k,folder in sources.items():
        checked_receipt(folder)
        digest=sha(folder/'receipt.json')
        if digest!=produced['source_receipts'][k]['sha256']:
            raise ValueError('Source receipt changed')
        pins[str(folder/'receipt.json')]=digest
    matrix=fasta(sources['matrix']/'matrix.faa')
    sites=table(sources['matrix']/'site_mapping.tsv')
    by_column={int(r['matrix_column_1based']):r for r in sites}
    columns=defaultdict(list)
    for row in sites:
        columns[row['marker']].append(int(row['matrix_column_1based']))
    if len(by_column)!=len(sites) or sorted(by_column)!=list(range(1,len(next(iter(matrix.values())))+1)):
        raise ValueError('Invalid source site map')
    encoded={}
    for row in table(sources['encodings']/'model_summary.tsv'):
        path=ROOT/row['encoding_path'];name=row['model_name']
        if name in encoded or sha(path)!=row['encoding_sha256']:
            raise ValueError('Duplicate or changed encoding')
        with np.load(path,allow_pickle=False) as d:
            seq=str(d['sequence']);states=str(d['states'])
            valid=d['valid'].astype(bool)
            confidence=d['feature_min_plddt'];error=d['feature_max_pae']
            accepted=valid & np.isfinite(confidence) & (confidence>=70) & np.isfinite(error) & (error<=10)
            ca=d['ca_plddt'].copy()
        if hashlib.sha256(seq.encode()).hexdigest()!=row['sequence_sha256'] or len(states)!=len(seq) or accepted.shape!=(len(seq),):
            raise ValueError('Encoding sequence or dimensions differ')
        encoded[name]=(seq,states,accepted,ca)
    links={}
    for row in table(sources['snapshot']/'marker_structure_links.tsv'):
        key=row['marker'],row['taxon_id']
        if key in links or key[0] not in columns or key[1] not in matrix:
            raise ValueError('Repeated or unknown model link')
        links[key]=Path(row['model_path']).stem
    observed=defaultdict(dict);seen=set()
    with gzip.open(sources['snapshot']/'matrix_to_structure_residues.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['marker'],row['taxon_id'];col=int(row['matrix_column_1based']);pos=int(row['protein_residue_1based'])-1
            identity=key+(col,)
            if identity in seen or key not in links or col not in by_column or by_column[col]['marker']!=key[0]:
                raise ValueError('Repeated or invalid residue mapping')
            seen.add(identity)
            seq,states,accepted,confidence=encoded[links[key]]
            if not 0<=pos<len(seq):
                raise ValueError('Residue out of bounds')
            aa=matrix[key[1]][col-1]
            if aa!=seq[pos] or not math.isclose(float(row['ca_plddt']),float(confidence[pos]),rel_tol=1e-5,abs_tol=1e-8):
                raise ValueError('Residue correspondence differs')
            if aa in ALPHABET and accepted[pos]:
                if states[pos] not in ALPHABET:
                    raise ValueError('Invalid observed structural symbol')
                observed[key][col]=(aa,states[pos])
    del seen,encoded
    coverage=table(a.inputs/'taxon_coverage.tsv');coverage_by_key={ (r['marker'],r['taxon_id']):r for r in coverage }
    if len(coverage_by_key)!=len(coverage) or set(coverage_by_key)!={(m,t) for m in columns for t in matrix}:
        raise ValueError('Incomplete marker/taxon coverage table')
    summary=table(a.inputs/'marker_summary.tsv');summary_by_marker={r['marker']:r for r in summary}
    if len(summary_by_marker)!=len(summary) or set(summary_by_marker)!=set(columns):
        raise ValueError('Marker universe differs')
    ready=set();taxa=set();cells=characters=retained_total=0
    for marker,cols in columns.items():
        threshold=max(50,math.ceil(len(cols)*0.3))
        eligible={t for t in matrix if len(observed.get((marker,t),{}))>=threshold}
        kept=sorted(set().union(*(set(observed[(marker,t)]) for t in eligible))) if eligible else []
        for t in matrix:
            row=coverage_by_key[marker,t]
            if (int(row['observed'])!=len(observed.get((marker,t),{}))
                    or row['taxon_eligible']!=str(t in eligible)
                    or int(row['required_observed_columns'])!=threshold):
                raise ValueError('Coverage or eligibility differs')
        row=summary_by_marker[marker]
        status='ready_for_inference' if len(eligible)>=4 else 'insufficient_structural_coverage'
        if int(row['eligible_taxa'])!=len(eligible) or int(row['retained_columns'])!=len(kept) or row['status']!=status:
            raise ValueError('Marker eligibility differs')
        folder=a.inputs/marker
        if len(eligible)<4:
            if folder.exists():raise ValueError('Ineligible marker emitted')
            continue
        ready.add(marker);taxa.update(eligible);cells+=len(eligible);retained_total+=len(kept)
        expected_columns=[dict(paired_column_1based=str(i+1),**by_column[c]) for i,c in enumerate(kept)]
        if table(folder/'columns.tsv')!=expected_columns:
            raise ValueError('Emitted column correspondence differs')
        for index,channel in enumerate(['aa','3di']):
            actual=fasta(folder/(channel+'.faa'))
            if set(actual)!=eligible:
                raise ValueError('Emitted taxon universe differs')
            for taxon,text in actual.items():
                expected=''.join(observed[(marker,taxon)].get(c,('?','?'))[index] for c in kept)
                if text!=expected:
                    raise ValueError('Emitted character or mask differs')
                if index==0:characters+=sum(c!='?' for c in text)
    if len(ready)!=produced['ready_markers'] or len(matrix)!=produced['taxa_audited'] or len(columns)!=produced['markers_audited']:
        raise ValueError('Receipt scope differs')
    for path,digest in pins.items():
        if sha(Path(path))!=digest:raise ValueError('Receipt changed during audit')
    checked_receipt(a.inputs)
    result=dict(status='passed_complete_paired_inputs_from_qualified_arrays_readback',markers=len(ready),taxa=len(taxa),marker_taxon_cells=cells,
                observed_paired_cells=characters,retained_marker_columns=retained_total,taxa_audited=len(matrix),markers_audited=len(columns),
                source_receipt_sha256=sha(a.inputs/'receipt.json'),script_sha256=sha(Path(__file__)),input_receipts=pins,
                scope='Every emitted character, observation mask, taxon eligibility and retained column reconstructed independently from mapped qualified arrays and source AA matrix; all 526-by-125 coverage cells checked. Does not re-infer structural features or establish biological accuracy, homology or branch rates.')
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
