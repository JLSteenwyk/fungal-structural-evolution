#!/usr/bin/env python3
"""Freeze lineage/length/confidence-stratified exact-sequence predictor controls."""
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable control input')
    mapping=ROOT/'results/structural_markers/gdm-expanded-v1'
    source=ROOT/'data/domains/marker-inputs-v1'
    checked_receipt(mapping);checked_receipt(source)
    manifest=ROOT/'metadata/analysis_manifest.tsv'
    taxa={r['taxon_id']:r for r in read_table(manifest)}
    models=json.loads((mapping/'model_provenance.json').read_text())
    by_model={(m['model_id'],str(m['version'])):m for m in models}
    sequences={hashlib.sha256(str(r.seq).encode()).hexdigest():str(r.seq) for r in SeqIO.parse(source/'sequences.faa','fasta')}
    links=read_table(mapping/'marker_structure_links.tsv')
    strata=defaultdict(set)
    for r in links:
        m=by_model[(r['model_id'],r['model_version'])]
        sequence=sequences.get(m['sequence_sha256'])
        if sequence is None:raise ValueError('Missing exact source sequence')
        if len(sequence)!=m['length']:raise ValueError('Sequence length differs')
        if len(sequence)>512 or not set(sequence)<=set('ACDEFGHIKLMNPQRSTVWY'):continue
        t=taxa[r['taxon_id']]
        group=(t['study_role'],t['lineage'].split(';')[0],(len(sequence)-1)//128+1,'below70' if m['mean_ca_plddt']<70 else 'at_least70')
        strata[group].add(m['sequence_sha256'])
    selected=set();selection=[]
    for group,pool in sorted(strata.items()):
        chosen=sorted(pool,key=lambda s:hashlib.sha256(('predictor-control-v1:'+s).encode()).hexdigest())[:4]
        selected.update(chosen)
        for sid in chosen:
            selection.append(dict(study_role=group[0],lineage_group=group[1],length_bin_128=group[2],af_mean_plddt_bin=group[3],eligible_unique_sequences=len(pool),sequence_id='S'+sid))
    if not selected:raise ValueError('No eligible controls')
    a.output.mkdir(parents=True)
    with (a.output/'candidates.faa').open('w') as f:
        for sid in sorted(selected):f.write('>S'+sid+'\n'+sequences[sid]+'\n')
    write_table(a.output/'selection_strata.tsv',selection)
    chosen_links=[r for r in links if r['sequence_sha256'] in selected]
    write_table(a.output/'reference_links.tsv',chosen_links)
    chosen_models=[m for m in models if m['sequence_sha256'] in selected]
    (a.output/'reference_models.json').write_text(json.dumps(chosen_models,indent=2)+'\n')
    emitted={r.id:str(r.seq) for r in SeqIO.parse(a.output/'candidates.faa','fasta')}
    if emitted!={'S'+s:sequences[s] for s in selected}:raise ValueError('Control FASTA readback differs')
    result={'status':'complete_prediction_queue_preparation','purpose':'Independent-predictor exact-sequence controls within the full project; not a pilot or experimental accuracy benchmark.',
            'source_mapping_receipt_sha256':sha(mapping/'receipt.json'),'sequence_source_receipt_sha256':sha(source/'receipt.json'),
            'manifest_sha256':sha(manifest),'script_sha256':sha(Path(__file__)),
            'prediction_candidates':len(selected),'candidate_residues':sum(len(sequences[s]) for s in selected),
            'candidate_taxa':len({r['taxon_id'] for r in chosen_links}),'reference_models':len(chosen_models),
            'reference_marker_links':len(chosen_links),'represented_lineage_role_groups':len({(r['study_role'],r['lineage_group']) for r in selection}),
            'selection':'Up to four unique sequences per role/manifest-lineage/128-residue-length/AF mean pLDDT (<70,>=70) stratum, smallest salted SHA256 ranks. Canonical full proteins <=512 residues only. Union sequences once; retain every associated marker/taxon/model identity. Same sequence may represent multiple strata, not independent replications.',
            'limitations':'Availability- and length-limited predictor sensitivity; no experimental reference. Both predictors use sequence, so agreement does not eliminate shared training or sequence-derived circularity. No ESMFold prediction or structural comparison completed by this stage.',
            'artifacts':{x.name:sha(x) for x in sorted(a.output.iterdir()) if x.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
