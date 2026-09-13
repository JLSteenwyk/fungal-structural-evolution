#!/usr/bin/env python3
"""Audit changes between frozen marker candidate queues without resubmitting old work."""
import argparse
import csv
import json
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from compare_marker_structures import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['previous','current','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    old_receipt=checked_receipt(args.previous); new_receipt=checked_receipt(args.current)
    if old_receipt['marker_input_receipt_sha256']!=new_receipt['marker_input_receipt_sha256']:
        raise ValueError('Different marker source universes')
    if args.output.exists():raise FileExistsError('Use a new immutable delta directory')
    old={r.id:str(r.seq) for r in SeqIO.parse(args.previous/'candidates.faa','fasta')}
    new={r.id:str(r.seq) for r in SeqIO.parse(args.current/'candidates.faa','fasta')}
    for key in old.keys()&new.keys():
        if old[key]!=new[key]:raise ValueError('Changed sequence under same identifier')
    states={r['sequence_id']:r['reuse_state'] for r in csv.DictReader((args.current/'all_marker_links.tsv').open(),delimiter='\t')}
    added=sorted(new.keys()-old.keys());removed=sorted(old.keys()-new.keys())
    args.output.mkdir(parents=True)
    with (args.output/'additional_candidates.faa').open('w') as out:
        for key in added:out.write('>'+key+'\n'+new[key]+'\n')
    rows=[]
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    for key in added:
        seq=new[key]; valid=not(set(seq)-canonical)
        rows.append({'sequence_id':key,'length':len(seq),'canonical_amino_acids':valid,
                     'within_existing_512_residue_limit':len(seq)<=512,'short_configuration_eligible':valid and len(seq)<=512})
    with (args.output/'additional_candidates.tsv').open('w') as out:
        writer=csv.DictWriter(out,['sequence_id','length','canonical_amino_acids','within_existing_512_residue_limit','short_configuration_eligible'],delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    with (args.output/'previous_candidates_now_reuse.tsv').open('w') as out:
        writer=csv.writer(out,delimiter='\t',lineterminator='\n');writer.writerow(['sequence_id','current_reuse_state'])
        for key in removed:
            if states[key] not in ['nominated_reuse_pending','verified_reuse_receipt']:raise ValueError('Unexplained removed candidate')
            writer.writerow([key,states[key]])
    receipt={'status':'complete_prediction_queue_delta','previous_receipt_sha256':sha(args.previous/'receipt.json'),
        'current_receipt_sha256':sha(args.current/'receipt.json'),'previous_candidates':len(old),'current_candidates':len(new),
        'additional_candidates':len(added),'additional_short_canonical_candidates':sum(r['short_configuration_eligible'] for r in rows),
        'previous_candidates_now_reuse':len(removed),'retained_previous_candidates':len(old.keys()&new.keys()),
        'script_sha256':sha(Path(__file__)),
        'interpretation':'This delta excludes all sequences already in the initial candidate queue, including its deferred long sequences. It is preparation, not a launched prediction run. Existing predictions and ongoing work must still be validated/reused when constructing the next execution queue; absence from the original queue is not proof that no model exists elsewhere.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
