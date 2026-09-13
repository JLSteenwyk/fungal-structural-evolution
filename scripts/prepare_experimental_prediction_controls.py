#!/usr/bin/env python3
"""Prepare all missing short ESMFold controls with exact experimental reference candidates."""
import argparse,hashlib,json
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT,sha,read_table
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use an immutable new queue')
    refs=ROOT/'results/experimental_structures/reference-metadata-v2';seqs=ROOT/'data/domains/marker-inputs-v1'
    checked_receipt(refs);checked_receipt(seqs)
    entries={r['entry_id']:r for r in read_table(refs/'entries.tsv')};targets=defaultdict(list)
    for r in read_table(refs/'entities.tsv'):
        if entries[r['entry_id']]['methodology']=='experimental':targets['S'+r['model_sequence_sha256']].append(r)
    sequences={r.id:str(r.seq) for r in SeqIO.parse(seqs/'sequences.faa','fasta')}
    reserved=defaultdict(list);pins={}
    for name in ['markers-v1','markers-followon-v1','predictor-controls-v1','ecology-markers-v1']:
        b=ROOT/'data/prediction_inputs'/name;checked_receipt(b);pins[name]=sha(b/'receipt.json')
        for r in SeqIO.parse(b/'candidates.faa','fasta'):reserved[r.id].append(name)
    disposition=[];links=[];chosen=set()
    for sid,hits in sorted(targets.items()):
        seq=sequences[sid]
        if sid!='S'+hashlib.sha256(seq.encode()).hexdigest():raise ValueError('Sequence hash differs')
        status=('deferred_noncanonical' if not set(seq)<=set('ACDEFGHIKLMNPQRSTVWY') else 'deferred_length' if len(seq)>512 else 'reserved_existing_queue' if sid in reserved else 'new_experimental_reference_control')
        if status=='new_experimental_reference_control':chosen.add(sid)
        disposition.append({'sequence_id':sid,'length':len(seq),'status':status,'existing_queues':';'.join(reserved.get(sid,[])),'reference_models':len({r['model_id'] for r in hits}),'reference_entries':len({r['entry_id'] for r in hits})})
        for h in hits:links.append({'sequence_id':sid,'predicted_model_id':h['model_id'],'entry_id':h['entry_id'],'entity_id':h['entity_id'],'sequence_sha256':h['model_sequence_sha256']})
    a.output.mkdir(parents=True)
    with (a.output/'candidates.faa').open('w') as f:
        for sid in sorted(chosen,key=lambda s:(len(sequences[s]),s)):f.write('>'+sid+'\n'+sequences[sid]+'\n')
    write_table(a.output/'sequence_disposition.tsv',disposition);write_table(a.output/'experimental_reference_links.tsv',links)
    emitted={r.id:str(r.seq) for r in SeqIO.parse(a.output/'candidates.faa','fasta')}
    if emitted!={s:sequences[s] for s in chosen} or chosen&reserved.keys():raise ValueError('Queue readback/disjointness failed')
    r={'status':'complete_prediction_queue_preparation','reference_receipt_sha256':sha(refs/'receipt.json'),'sequence_source_receipt_sha256':sha(seqs/'receipt.json'),'prior_input_receipt_sha256':pins,'script_sha256':sha(Path(__file__)),'experimental_reference_sequences':len(targets),'prediction_candidates':len(chosen),'candidate_residues':sum(len(sequences[s]) for s in chosen),'sequence_status_counts':dict(Counter(r['status'] for r in disposition)),'selection':'All exact full canonical experimental-reference sequences eligible for the existing <=512-residue ESMFold configuration, excluding previously reserved sequence IDs. No selection on structure agreement, experimental resolution, or phenotype. No truncation.','limitations':'Preparation only; prior reservation is not successful prediction. Existing hardware is busy, no additional GPU process launched. Three-way AF/ESMFold/experiment comparison, experimental quality/context and training/template independence remain pending; reference sampling is concentrated in Ascomycota.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
