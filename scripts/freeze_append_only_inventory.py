#!/usr/bin/env python3
"""Freeze a verified complete-line prefix of a concurrently appended JSONL log."""
import argparse
import hashlib
import json
from pathlib import Path


def prefix_hash(path, size):
    h=hashlib.sha256()
    with path.open('rb') as f:
        remaining=size
        while remaining:
            block=f.read(min(8388608,remaining))
            if not block:raise ValueError('Source truncated')
            h.update(block);remaining-=len(block)
    return h.hexdigest()


def freeze(source, output):
    output.mkdir(parents=True,exist_ok=False)
    before=source.stat();target=output/'inventory.jsonl'
    remaining=before.st_size;buffer=b'';lines=0;written=0;digest=hashlib.sha256()
    with source.open('rb') as src,target.open('wb') as dest:
        while remaining:
            block=src.read(min(8388608,remaining))
            if not block:raise ValueError('Source truncated during freeze')
            remaining-=len(block);buffer+=block
            boundary=buffer.rfind(b'\n')+1
            if boundary:
                complete,buffer=buffer[:boundary],buffer[boundary:]
                dest.write(complete);digest.update(complete);written+=len(complete);lines+=complete.count(b'\n')
    after=source.stat()
    if before.st_ino!=after.st_ino or before.st_dev!=after.st_dev or after.st_size<before.st_size:
        raise ValueError('Source identity changed or was truncated')
    if not lines or prefix_hash(source,written)!=digest.hexdigest() or prefix_hash(target,written)!=digest.hexdigest():
        raise ValueError('Frozen prefix failed source readback')
    receipt=dict(status='complete_verified_append_only_inventory_prefix',source=str(source),source_initial_bytes=before.st_size,
                 source_observed_final_bytes=after.st_size,frozen_bytes=written,complete_lines=lines,
                 incomplete_tail_bytes_excluded=len(buffer),inventory_sha256=digest.hexdigest(),
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 scope='Exact complete-line prefix at observed initial byte length, verified by rereading the source prefix. JSON content and model validity require downstream checks. Concurrent later appends are excluded.')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();print(json.dumps(freeze(a.source,a.output),indent=2))


if __name__=='__main__':main()
