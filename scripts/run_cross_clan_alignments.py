#!/usr/bin/env python3
"""Checkpoint both input orders of every cross-clan candidate interval alignment."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time
from Bio import SeqIO


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def parse_output(raw, sequences):
    lengths = [int(re.search(r'Length of Structure_'+str(i)+r':\s+(\d+) residues', raw).group(1)) for i in [1,2]]
    match = re.search(r'Aligned length=\s*(\d+), RMSD=\s*([0-9.]+), Seq_ID=n_identical/n_aligned=\s*([0-9.]+)', raw)
    if match is None: raise ValueError('Missing native summary')
    n, rmsd, identity = int(match[1]), float(match[2]), float(match[3])
    scores = [float(re.search(r'TM-score=\s*([0-9.]+) \(normalized by length of Structure_'+str(i)+':', raw).group(1)) for i in [1,2]]
    lines = raw.splitlines(); index = next(i for i,l in enumerate(lines) if l.startswith('(":" denotes residue pairs'))
    left, marks, right = lines[index+1:index+4]
    if len(left) != len(right) or len(marks) != len(left): raise ValueError('Alignment widths differ')
    if [left.replace('-',''), right.replace('-','')] != sequences: raise ValueError('Native alignment sequences differ')
    if lengths != list(map(len,sequences)): raise ValueError('Native lengths differ')
    if not 0 < n <= min(lengths) or not 0 <= identity <= 1 or rmsd < 0: raise ValueError('Invalid alignment summary')
    if not all(math.isfinite(v) for v in [rmsd,identity,*scores]) or not all(0 <= v <= 1 for v in scores): raise ValueError('Invalid score')
    return dict(length_left=lengths[0],length_right=lengths[1],aligned_length=n,rmsd=rmsd,
                sequence_identity=identity,tm_left=scores[0],tm_right=scores[1],
                alignment_left=left,alignment_marks=marks,alignment_right=right)


def main():
    ap = argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Plan changed')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Input changed: '+p)
    verify()
    coords=Path(plan['coordinates']); cr=json.loads((coords/'receipt.json').read_text())
    conf=json.loads(Path(plan['confidence_receipt']).read_text())
    candidate=json.loads(Path(plan['candidate_receipt']).read_text())
    if cr['status']!='complete_cross_clan_coordinate_materialization' or conf['status']!='complete_cross_clan_confidence_manifest_join' or candidate['status']!='complete_cross_clan_candidate_inventory':raise ValueError('Incomplete input')
    if candidate['artifacts']['candidate_members.tsv']!=plan['pins'][plan['members']]:raise ValueError('Candidate binding differs')
    if conf['source_hashes'][plan['members']]!=plan['pins'][plan['members']]:raise ValueError('Confidence candidate binding differs')
    for name,h in cr['pdb_hashes'].items():
        if sha(coords/'pdb'/name)!=h:raise ValueError('Coordinate changed')
    if sha(coords/'sequences.faa')!=cr['artifacts']['sequences.faa']:raise ValueError('Sequence artifact changed')
    sequences={s.id:str(s.seq) for s in SeqIO.parse(coords/'sequences.faa','fasta')}
    groups=defaultdict(lambda:[set(),set()])
    with open(plan['members']) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['model_exclusive_within_pair_cluster']=='1':
                key=tuple(r[k] for k in ['boundary','representative','pfam_left','pfam_right'])
                groups[key][int(r['side']=='right')].add(r['interval_id'])
    pairs={tuple(sorted((a,b))) for left,right in groups.values() for a in left for b in right}
    if len(pairs)!=conf['workloads_by_fraction_ge70']['0']['unique_unordered_interval_pairs']:raise ValueError('Pair scope differs')
    jobs=[pair for a,b in sorted(pairs) for pair in [(a,b),(b,a)]]
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    (out/'pairs').mkdir(exist_ok=True)
    def run(pair):
        a,b=pair;key=hashlib.sha256((a+'\t'+b).encode()).hexdigest();dest=out/'pairs'/(key+'.json')
        command=[plan['usalign'],str(coords/'pdb'/(a+'.pdb')),str(coords/'pdb'/(b+'.pdb')),*plan['options']]
        if dest.exists():
            r=json.loads(dest.read_text())
            if r['plan_sha256']!=ph or r['intervals']!=[a,b] or r['command']!=command:raise ValueError('Checkpoint binding differs')
            if parse_output(r['stdout'],[sequences[a],sequences[b]])!=r['metrics']:raise ValueError('Checkpoint metrics differ')
            return
        started=time.monotonic();result=subprocess.run(command,capture_output=True,text=True,check=True,timeout=plan['per_pair_timeout_seconds'])
        metrics=parse_output(result.stdout,[sequences[a],sequences[b]])
        record=dict(plan_sha256=ph,intervals=[a,b],command=command,metrics=metrics,
                    stdout=result.stdout,stderr=result.stderr,elapsed_seconds=time.monotonic()-started)
        temp=dest.with_suffix('.tmp');temp.write_text(json.dumps(record,indent=2)+'\n');temp.replace(dest)
    with ThreadPoolExecutor(max_workers=plan['resources']['cpus']) as pool:
        for n,_ in enumerate(pool.map(run,jobs),1):
            if n%100==0:print(n,'/',len(jobs),'directed alignments completed',flush=True)
    verify()
    for name,h in cr['pdb_hashes'].items():
        if sha(coords/'pdb'/name)!=h:raise ValueError('Coordinate changed during run')
    files=sorted((out/'pairs').glob('*.json'))
    if len(files)!=len(jobs):raise ValueError('Unexpected checkpoint count')
    receipt=dict(status='complete_cross_clan_alignments_pending_independent_readback',plan_sha256=ph,
                 unordered_pairs=len(pairs),directed_alignments=len(jobs),
                 artifacts={str(p.relative_to(out)):sha(p) for p in files},
                 scope='Full candidate interval grid, both input orders, monomeric protein structural alignment. '
                 'Raw output and alignment sequences retained. Confidence stratification and independent numeric readback pending. '
                 'No PAE mask, phylogenetic test, established homology or functional inference.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(receipt['status'],flush=True)


if __name__=='__main__':main()
