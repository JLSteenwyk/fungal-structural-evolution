#!/usr/bin/env python3
"""Partition a cleanly stopped ESMFold queue without changing completed models."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import psutil


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','predictions','checkpoint_record','output']:
        ap.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    a=ap.parse_args()
    stopped=json.loads(a.checkpoint_record.read_text())
    try:
        process=psutil.Process(stopped['pid'])
        if process.create_time()==stopped['created'] and process.status()!=psutil.STATUS_ZOMBIE:
            raise RuntimeError('Original predictor still alive')
    except psutil.NoSuchProcess:
        pass
    chunk=json.loads((a.predictions/'last_chunk.json').read_text())
    config=json.loads((a.predictions/'config.json').read_text())
    config_sha=sha(a.predictions/'config.json')
    if not chunk['interrupted'] or chunk['oom_deferred'] or chunk['config_sha256']!=config_sha:
        raise ValueError('Expected clean interrupted chunk without OOM')
    receipt=json.loads((a.inputs/'receipt.json').read_text())
    if sha(a.inputs/'receipt.json')!=config['input_receipt_sha256']:
        raise ValueError('Input receipt mismatch')
    for name,digest in receipt['artifacts'].items():
        if sha(a.inputs/name)!=digest:
            raise ValueError('Changed input artifact')
    sequences={}; sid=None; parts=[]
    for line in (a.inputs/'candidates.faa').read_text().splitlines():
        if line.startswith('>'):
            if sid is not None:
                if sid in sequences:raise ValueError('Duplicate input ID')
                sequences[sid]=''.join(parts)
            sid=line[1:].strip();parts=[]
        else:parts.append(line.strip())
    if sid is not None:
        if sid in sequences:raise ValueError('Duplicate input ID')
        sequences[sid]=''.join(parts)
    for sid,sequence in sequences.items():
        if sid!='S'+hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError('Sequence digest mismatch')
    complete={}
    for path in sorted(a.predictions.glob('S*.json')):
        row=json.loads(path.read_text());sid=row['sequence_id']
        if row['status']!='verified_prediction' or sid not in sequences or row['config_sha256']!=config_sha:
            raise ValueError('Invalid completed model receipt')
        if row['length']!=len(sequences[sid]) or row['sequence_sha256']!=sid[1:]:
            raise ValueError('Completed sequence mismatch')
        for name,digest in row['artifacts'].items():
            if sha(a.predictions/name)!=digest:
                raise ValueError('Changed completed prediction')
        complete[sid]=sha(path)
    remaining=set(sequences)-set(complete)
    if len(remaining)!=chunk['remaining_eligible'] or len(complete)!=chunk['new_predictions']+chunk['cached_predictions']:
        raise ValueError('Checkpoint count mismatch')
    queues=[[],[]];cost=[0,0]
    for sid in sorted(remaining,key=lambda s:(-len(sequences[s]),s)):
        index=0 if cost[0]<=cost[1] else 1
        queues[index].append(sid);cost[index]+=len(sequences[sid])**3
    if set(queues[0])&set(queues[1]) or set(queues[0])|set(queues[1])!=remaining:
        raise ValueError('Split coverage or overlap error')
    with (a.inputs/'all_marker_links.tsv').open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;links=list(reader)
    a.output.mkdir(parents=True,exist_ok=False)
    reports=[]
    for index,queue in enumerate(queues):
        folder=a.output/f'gpu{index}';folder.mkdir()
        (folder/'candidates.faa').write_text(''.join('>'+s+'\n'+sequences[s]+'\n' for s in queue))
        selected=set(queue);selected_links=[r for r in links if r['sequence_id'] in selected]
        if {r['sequence_id'] for r in selected_links}!=selected:
            raise ValueError('Missing sequence links')
        with (folder/'all_marker_links.tsv').open('w') as handle:
            w=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected_links)
        r=dict(status='complete_disjoint_remaining_queue',gpu_index=index,prediction_candidates=len(queue),
               marker_links=len(selected_links),residues=sum(len(sequences[s]) for s in queue),
               cubic_length_cost=cost[index],original_input_receipt_sha256=sha(a.inputs/'receipt.json'),
               original_chunk_sha256=sha(a.predictions/'last_chunk.json'),
               artifacts={name:sha(folder/name) for name in ['candidates.faa','all_marker_links.tsv']})
        (folder/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');reports.append(r)
    (a.output/'completed_source_receipts.json').write_text(json.dumps(complete,sort_keys=True)+'\n')
    result=dict(status='complete_two_gpu_remaining_partition',total=len(sequences),completed=len(complete),remaining=len(remaining),
                queues=reports,original_predictions=str(a.predictions.resolve()),script_sha256=sha(Path(__file__)),
                completed_source_receipts_sha256=sha(a.output/'completed_source_receipts.json'),
                interpretation='Completed original artifacts preserved and rehashed. Remaining full sequences divided into disjoint queues, balanced by cubic length; separate output provenance is required. Original wrapper is superseded, not a full completed prediction tier.')
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
