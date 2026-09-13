#!/usr/bin/env python3
"""Freeze missing same-method marker inputs for every curated ecological species."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable input directory')
    source=ROOT/'data/domains/marker-inputs-v1';checked_receipt(source)
    evidence=ROOT/'metadata/species_ecology_evidence.tsv'
    er=ROOT/'metadata/species_ecology_evidence_receipt.json'
    if sha(evidence)!=json.loads(er.read_text())['table_sha256']:raise ValueError('Changed ecology evidence')
    taxa={r['taxon_id'] for r in read_table(evidence)}
    sequences={r.id:str(r.seq) for r in SeqIO.parse(source/'sequences.faa','fasta')}
    links=[r for r in read_table(source/'protein_links.tsv') if r['taxon_id'] in taxa]
    # Existing immutable input queues reserve sequences whether completed or still pending.
    reserved=set();prior={}
    for name in ['markers-v1','markers-followon-v1','predictor-controls-v1']:
        b=ROOT/'data/prediction_inputs'/name;checked_receipt(b)
        prior[name]=sha(b/'receipt.json')
        reserved.update(r.id for r in SeqIO.parse(b/'candidates.faa','fasta'))
    rows=[];chosen=set()
    for sid in sorted({r['sequence_id'] for r in links}):
        sequence=sequences[sid]
        if sid!='S'+hashlib.sha256(sequence.encode()).hexdigest():raise ValueError('Sequence identity differs')
        status=('deferred_noncanonical' if not set(sequence)<=set('ACDEFGHIKLMNPQRSTVWY') else
                'deferred_length' if len(sequence)>512 else
                'reserved_existing_queue' if sid in reserved else 'new_same_method_candidate')
        rows.append({'sequence_id':sid,'length':len(sequence),'status':status})
        if status=='new_same_method_candidate':chosen.add(sid)
    a.output.mkdir(parents=True)
    with (a.output/'candidates.faa').open('w') as f:
        for sid in sorted(chosen,key=lambda s:(len(sequences[s]),s)):f.write('>'+sid+'\n'+sequences[sid]+'\n')
    write_table(a.output/'sequence_disposition.tsv',rows)
    write_table(a.output/'ecology_marker_links.tsv',links)
    emitted={r.id:str(r.seq) for r in SeqIO.parse(a.output/'candidates.faa','fasta')}
    if emitted!={s:sequences[s] for s in chosen} or chosen&reserved:raise ValueError('Queue readback/overlap failure')
    r={'status':'complete_prediction_queue_preparation','purpose':'Same-method ecological marker coverage, regardless of existing AlphaFold availability; all curated species, not phenotype-selected families.',
       'ecology_evidence_receipt_sha256':sha(er),'ecology_table_sha256':sha(evidence),'source_receipt_sha256':sha(source/'receipt.json'),
       'prior_input_receipt_sha256':prior,'script_sha256':sha(Path(__file__)),'prediction_candidates':len(chosen),
       'candidate_residues':sum(len(sequences[s]) for s in chosen),'curated_taxa':len(taxa),'candidate_taxa':len({r['taxon_id'] for r in links if r['sequence_id'] in chosen}),
       'sequence_status_counts':dict(Counter(r['status'] for r in rows)),
       'limitations':'Prepared only, not launched. Existing-queue reservation does not assert successful prediction. Full proteins only, no truncation; length/alphabet exclusions retained. Confidence, orthology and transition-replication gates still apply.',
       'artifacts':{x.name:sha(x) for x in sorted(a.output.iterdir()) if x.is_file()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
