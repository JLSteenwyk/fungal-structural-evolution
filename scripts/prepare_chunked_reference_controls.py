#!/usr/bin/env python3
"""Prepare the last unreserved reference and an explicit chunk-size comparison control."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from audit_joint_path_uncertainty import checked,rows,sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    source=Path('data/prediction_inputs/experimental-controls-v1');checked(source)
    prior=Path('results/predictions/experimental-controls-513-768-audit-v1');checked(prior)
    control=max(rows(prior/'predictions.tsv'),key=lambda r:int(r['length']))
    sequences={r.id:str(r.seq) for r in SeqIO.parse('data/domains/marker-inputs-v1/sequences.faa','fasta')}
    reserved=set();pins={}
    for f in Path('data/prediction_inputs').glob('*/candidates.faa'):
        r=json.loads((f.parent/'receipt.json').read_text())
        if sha(f)!=r['artifacts']['candidates.faa']:raise ValueError('Changed input queue')
        pins[str(f)]=sha(f);reserved.update(r.id for r in SeqIO.parse(f,'fasta'))
    remaining={r['sequence_id'] for r in rows(source/'sequence_disposition.tsv') if r['sequence_id'] not in reserved}
    if len(remaining)!=1 or len(sequences[next(iter(remaining))])!=2413:raise ValueError('Remaining reference universe changed')
    chosen=remaining|{control['sequence_id']}
    if len(chosen)!=2:raise ValueError('Invalid comparison control')
    for sid in chosen:
        if sid!='S'+hashlib.sha256(sequences[sid].encode()).hexdigest() or not set(sequences[sid])<=set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('Invalid full canonical sequence')
    a.output.mkdir(parents=True)
    (a.output/'candidates.faa').write_text(''.join('>'+sid+'\n'+sequences[sid]+'\n' for sid in sorted(chosen,key=lambda s:len(sequences[s]))))
    links=[r for r in rows(source/'experimental_reference_links.tsv') if r['sequence_id'] in chosen]
    if {r['sequence_id'] for r in links}!=chosen:raise ValueError('Missing reference links')
    write_table(a.output/'experimental_reference_links.tsv',links)
    write_table(a.output/'sequence_disposition.tsv',[dict(sequence_id=sid,length=len(sequences[sid]),role='new_long_reference' if sid in remaining else 'intentional_chunk_size_sensitivity_control') for sid in sorted(chosen)])
    result=dict(status='complete_prediction_queue_preparation',prediction_candidates=2,new_reference_sequences=1,intentional_repredictions=1,
        original_reference_receipt_sha256=sha(source/'receipt.json'),comparison_control_audit_sha256=sha(prior/'receipt.json'),
        comparison_control_prediction_receipt_sha256=control['prediction_receipt_sha256'],comparison_control_sequence_id=control['sequence_id'],
        reservation_fasta_sha256=pins,script_sha256=sha(Path(__file__)),
        interpretation='New 2413-residue full reference plus longest successfully audited 513–768 control, selected by length rather than agreement. The existing control is deliberately repeated under chunk16 to assess numerical sensitivity; not pooled with original predictions.',
        artifacts={f.name:sha(f) for f in a.output.iterdir()})
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print('Prepared full lengths',sorted(len(sequences[s]) for s in chosen))


if __name__=='__main__':main()
