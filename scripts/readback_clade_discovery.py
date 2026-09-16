#!/usr/bin/env python3
"""Independently validate exported clade partitions against native clusters and FASTAs."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from orthofinder.tools import mcl


def sha(p):
    with Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--synthetic-exact-copies',action='store_true')
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    plan=json.loads(a.plan.read_text());root=Path(plan['output']);overall=json.loads((root/'receipt.json').read_text())
    if overall['plan_sha256']!=sha(a.plan) or overall['clades']!=len(plan['cases']):raise ValueError('Overall receipt mismatch')
    reports=[]
    for case in plan['cases']:
        key=case['clade_id'];folder=root/key;r=json.loads((folder/'receipt.json').read_text())
        if sha(folder/'receipt.json')!=overall['clade_receipts'][key] or r['plan_sha256']!=sha(a.plan):raise ValueError('Changed clade receipt')
        for name,h in r['artifacts'].items():
            if sha(folder/name)!=h:raise ValueError('Changed output '+name)
        expected=set();copy_groups=defaultdict(set)
        for item in case['input_files']:
            path=Path(case['input_directory'])/item['filename']
            if sha(path)!=item['sha256']:raise ValueError('Input changed')
            for rec in SeqIO.parse(path,'fasta'):
                if rec.id in expected:raise ValueError('Duplicate source identity')
                expected.add(rec.id)
                if a.synthetic_exact_copies:copy_groups[str(rec.seq)].add(rec.id)
        clusters=[folder/n for n in r['artifacts'] if n.endswith('_id_pairs.txt')]
        idsfiles=[folder/n for n in r['artifacts'] if n.endswith('SequenceIDs.txt')]
        if len(clusters)!=1 or len(idsfiles)!=1:raise ValueError('Ambiguous artifacts')
        ids=dict(line.split(': ',1) for line in idsfiles[0].read_text().splitlines())
        native=mcl.GetPredictedOGs(str(clusters[0]));ngroups=[{ids[x] for x in g} for g in native]
        exported=defaultdict(set);seen=set()
        with (folder/'family_membership.tsv').open() as h:
            for row in csv.DictReader(h,delimiter='\t'):
                gene=row['original_native_id']
                if gene in seen:raise ValueError('Repeated export gene')
                seen.add(gene);exported[row['local_family']].add(gene)
        if set(exported)!={f'OG{i:07d}' for i in range(len(ngroups))}:raise ValueError('Family indexing differs')
        if any(exported[f'OG{i:07d}']!=g for i,g in enumerate(ngroups)):raise ValueError('Native/export partition differs')
        if seen!=expected or len(seen)!=r['proteins'] or len(exported)!=r['families']:raise ValueError('Protein coverage differs')
        if a.synthetic_exact_copies and {frozenset(g) for g in ngroups}!={frozenset(g) for g in copy_groups.values()}:raise ValueError('Synthetic exact-copy classes differ')
        reports.append(dict(clade_id=key,proteins=len(seen),families=len(exported),synthetic_exact_copies_checked=a.synthetic_exact_copies))
    result=dict(status='passed_complete_native_clade_discovery_readback',clades=reports,plan_sha256=sha(a.plan),receipt_sha256=sha(root/'receipt.json'),script_sha256=sha(Path(__file__)),native_reader_sha256=sha(Path(mcl.__file__)),scope='Independent BioPython FASTA parsing and native cluster reader check every exported family and original protein ID. Does not establish homology accuracy, cross-clade completeness or reconciled orthology.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
