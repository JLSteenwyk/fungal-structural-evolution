#!/usr/bin/env python3
"""Materialize every cross-clan candidate domain with full member-byte checks."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import tarfile
from Bio.Data.PDBData import protein_letters_3to1


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def validate(blob, row):
    if hashlib.sha256(blob).hexdigest() != row['pdb_sha256']: raise ValueError('Member checksum differs')
    ca = {}; atoms = set()
    for line in blob.decode().splitlines():
        if not line.startswith('ATOM  '):
            if line not in ['TER', 'END']: raise ValueError('Unsupported record')
            continue
        pos = int(line[22:26]); atom = line[12:16].strip()
        if (pos, atom) in atoms or line[21] != 'A': raise ValueError('Duplicate atom or unexpected chain')
        atoms.add((pos, atom))
        values = [float(line[a:b]) for a,b in [(30,38),(38,46),(46,54),(60,66)]]
        if not all(math.isfinite(x) for x in values) or not 0 <= values[-1] <= 100: raise ValueError('Invalid coordinate/confidence')
        if atom == 'CA': ca[pos] = (protein_letters_3to1[line[17:20]], values[-1])
    n = int(row['residues'])
    if set(ca) != set(range(1,n+1)): raise ValueError('Incomplete C-alpha grid')
    sequence = ''.join(ca[i][0] for i in range(1,n+1))
    if hashlib.sha256(sequence.encode()).hexdigest() != row['fragment_sequence_sha256']: raise ValueError('Fragment sequence differs')
    mean = sum(ca[i][1] for i in ca)/n
    if abs(mean-float(row['mean_ca_plddt'])) > .005001: raise ValueError('Mean confidence differs beyond rounding')
    return sequence


def main():
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args(); plan = json.loads(args.plan.read_text()); ph = sha(args.plan)
    def verify():
        if sha(args.plan) != ph: raise ValueError('Plan changed')
        for p,d in plan['pins'].items():
            if sha(p) != d: raise ValueError('Source changed: '+p)
    verify()
    receipt = json.loads(Path(plan['confidence_receipt']).read_text())
    if receipt['status'] != 'complete_cross_clan_confidence_manifest_join': raise ValueError('Completed confidence join required')
    if receipt['artifacts']['interval_confidence.tsv'] != plan['pins'][plan['intervals']]: raise ValueError('Confidence table binding differs')
    with open(plan['intervals']) as f: rows = list(csv.DictReader(f, delimiter='\t'))
    ids = [r['interval_id'] for r in rows]
    if len(set(ids)) != len(ids) or len(ids) != receipt['candidate_intervals']: raise ValueError('Incomplete candidate grid')
    out = Path(plan['output']); out.mkdir(parents=True, exist_ok=False); (out/'pdb').mkdir()
    grouped = defaultdict(list)
    for row in rows:
        if not re.fullmatch(r'D[0-9a-f]{64}', row['interval_id']) or row['member'] != row['interval_id']+'.pdb': raise ValueError('Invalid interval identity')
        grouped[row['archive_path']].append(row)
    byte_count = 0
    with (out/'sequences.faa').open('w') as fasta:
        for archive_path, entries in sorted(grouped.items()):
            if shutil.disk_usage(out).free < plan['resources']['minimum_free_disk_gib']*2**30: raise ValueError('Disk reserve reached')
            with tarfile.open(archive_path, 'r:') as archive:
                wanted = {r['member']:r for r in entries}; seen = set()
                for member in archive:
                    if member.name not in wanted: continue
                    if member.name in seen or not member.isfile(): raise ValueError('Invalid candidate archive member')
                    seen.add(member.name); row = wanted[member.name]
                    blob = archive.extractfile(member).read(); sequence = validate(blob, row)
                    target = out/'pdb'/member.name; target.write_bytes(blob)
                    if sha(target) != row['pdb_sha256']: raise ValueError('Written member differs')
                    fasta.write('>'+row['interval_id']+'\n'+sequence+'\n'); byte_count += len(blob)
                if seen != set(wanted): raise ValueError('Missing archive member')
    verify()
    result = dict(status='complete_cross_clan_coordinate_materialization', plan_sha256=ph,
                  intervals=len(rows), archives=len(grouped), pdb_bytes=byte_count,
                  artifacts={'sequences.faa':sha(out/'sequences.faa')},
                  pdb_hashes={r['interval_id']+'.pdb':r['pdb_sha256'] for r in rows},
                  scope='Every candidate PDB member checked against audited manifest bytes, decoded for sequence, '
                  'C-alpha coverage, finite coordinates, and mean confidence within rounding tolerance. '
                  'Output bytes rehashed. Full archive bytes not rehashed; selected members checked individually. '
                  'No structural alignments or homology conclusions.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='pdb_hashes'},indent=2))


if __name__ == '__main__': main()
