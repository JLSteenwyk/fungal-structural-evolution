#!/usr/bin/env python3
"""Prepare a disjoint follow-on marker queue, checking newly available models."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import ROOT, sha, read_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'current', 'delta', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--existing-predictions', type=Path, nargs='+', required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('Use a new immutable input snapshot')
    for path in [args.previous, args.current, args.delta]: checked_receipt(path)
    delta = json.loads((args.delta/'receipt.json').read_text())
    if delta['previous_receipt_sha256'] != sha(args.previous/'receipt.json') or delta['current_receipt_sha256'] != sha(args.current/'receipt.json'):
        raise ValueError('Queue delta has different input provenance')
    previous = {r.id for r in SeqIO.parse(args.previous/'candidates.faa','fasta')}
    new_records = list(SeqIO.parse(args.delta/'additional_candidates.faa','fasta'))
    new = {r.id:str(r.seq) for r in new_records}
    if len(new) != len(new_records) or previous & new.keys(): raise ValueError('Repeated or previously queued sequence')
    for sid, sequence in new.items():
        if sid != 'S'+hashlib.sha256(sequence.encode()).hexdigest(): raise ValueError('Sequence identifier mismatch')
    raw = (ROOT/'data/raw/afdb_models.jsonl').read_bytes().splitlines(keepends=True)
    if raw and not raw[-1].endswith(b'\n'): raw.pop()
    frozen = b''.join(raw)
    latest = {}
    for line in frozen.decode().splitlines():
        row = json.loads(line); latest[row['uniprot_accession']] = row
    available = defaultdict(list)
    for row in latest.values():
        if row['status'] == 'verified':
            for model in row['models']:
                sid = 'S'+model['sequence_sha256']
                if sid in new:
                    if sha(ROOT/model['path']) != model['sha256']: raise ValueError('Changed reusable model')
                    available[sid].append({'path':model['path'],'sha256':model['sha256'], 'source':model.get('provider'), 'tool':model.get('tool')})
    for folder in args.existing_predictions:
        for sid in new:
            path = folder/(sid+'.json')
            if not path.exists(): continue
            row = json.loads(path.read_text())
            if row['status'] != 'verified_prediction' or row['sequence_sha256'] != sid[1:] or row['sequence_id'] != sid:
                raise ValueError('Invalid existing prediction')
            for name, checksum in row['artifacts'].items():
                if sha(folder/name) != checksum: raise ValueError('Changed existing prediction artifact')
            available[sid].append({'path':str(path), 'sha256':sha(path), 'source':'existing_local_prediction', 'tool':row['source']})
    links = [r for r in read_table(args.current/'all_marker_links.tsv') if r['sequence_id'] in new]
    if {r['sequence_id'] for r in links} != set(new): raise ValueError('Incomplete follow-on taxon linkage')
    dispositions = []
    canonical = set('ACDEFGHIKLMNPQRSTVWY')
    for sid, seq in sorted(new.items()):
        state = ('verified_model_available' if available[sid] else 'noncanonical_deferred' if not set(seq)<=canonical else 'length_deferred' if len(seq)>512 else 'short_prediction_eligible')
        dispositions.append({'sequence_id':sid, 'length':len(seq), 'queue_disposition':state,
            'existing_models_json':json.dumps(available[sid],sort_keys=True)})
    candidates = [r for r in dispositions if r['queue_disposition']!='verified_model_available']
    args.output.mkdir(parents=True)
    (args.output/'afdb_snapshot.jsonl').write_bytes(frozen)
    with (args.output/'candidates.faa').open('w') as handle:
        for row in candidates: handle.write('>'+row['sequence_id']+'\n'+new[row['sequence_id']]+'\n')
    for name, rows in [('all_marker_links.tsv',links),('queue_disposition.tsv',dispositions)]:
        with (args.output/name).open('w') as handle:
            writer=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    from collections import Counter
    result={'status':'complete_disjoint_followon_prediction_inputs','additional_sequences_screened':len(new),
        'prediction_candidates':len(candidates),'taxon_marker_links':len(links),'linked_taxa':len({r['taxon_id'] for r in links}),
        'dispositions':dict(Counter(r['queue_disposition'] for r in dispositions)),
        'previous_receipt_sha256':sha(args.previous/'receipt.json'),'current_receipt_sha256':sha(args.current/'receipt.json'),
        'delta_receipt_sha256':sha(args.delta/'receipt.json'),'script_sha256':sha(Path(__file__)),
        'existing_prediction_directories':[str(p) for p in args.existing_predictions],
        'interpretation':'Excludes every original candidate, including long/deferred ones; all follow-on taxon/protein links preserved. Current verified exact-sequence models and existing local predictions checked before preparation. Long/noncanonical sequences remain explicit. Input preparation does not launch predictions; recheck reuse and GPU availability before execution. Source inventory states in all_marker_links are historical, whereas queue_disposition records this snapshot.',
        'artifacts':{p.name:sha(p) for p in args.output.iterdir()}}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__': main()
