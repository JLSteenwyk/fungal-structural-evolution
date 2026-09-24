#!/usr/bin/env python3
"""Independently rebuild domain cluster joins and counts from full source rows."""
import argparse
import csv
import hashlib
import json
import sqlite3
import time
from collections import Counter, defaultdict
from pathlib import Path
import psutil


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def connect(path):
    return sqlite3.connect('file:'+str(Path(path).resolve())+'?mode=ro', uri=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    ph = sha(args.plan)

    def verify():
        if sha(args.plan) != ph:
            raise ValueError('Audit plan changed')
        for name, digest in plan['pins'].items():
            if sha(name) != digest:
                raise ValueError('Pinned input changed: '+name)

    verify()
    pred = plan['predecessor']
    while True:
        try:
            proc = psutil.Process(pred['pid'])
            if proc.create_time() != pred['create_time'] or proc.status() == psutil.STATUS_ZOMBIE:
                break
            if proc.cmdline() != pred['command']:
                raise ValueError('Producer command changed')
        except psutil.NoSuchProcess:
            break
        time.sleep(20)
    verify()
    p = json.loads(Path(plan['producer_plan']).read_text())
    root = Path(p['output'])
    rp = root/'receipt.json'
    r = json.loads(rp.read_text())
    receipt_hash = sha(rp)
    if r['status'] != 'complete_domain_cluster_composition_pending_independent_readback' or r['plan_sha256'] != sha(plan['producer_plan']):
        raise ValueError('Producer incomplete or mismatched')

    def verify_outputs():
        if sha(rp) != receipt_hash:
            raise ValueError('Receipt changed')
        for name, digest in r['artifacts'].items():
            if sha(root/name) != digest:
                raise ValueError('Output changed: '+name)

    verify_outputs()
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    db = connect(root/'domain_cluster_composition.sqlite')
    source = connect(p['source_database'])
    intervals = {}
    with Path(p['intervals']).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = row['interval_id']
            if key in intervals:
                raise ValueError('Duplicate interval')
            intervals[key] = (row['model_key'], int(row['start']), int(row['end']))
    remaining = dict(intervals)
    for key, model, start, end in db.execute('SELECT interval_id,model,start,end FROM intervals'):
        if remaining.pop(key, None) != (model,start,end):
            raise ValueError('Interval source differs')
    if remaining:
        raise ValueError('Missing interval')
    members = {}
    with Path(p['members']).open() as f:
        for rep, key in csv.reader(f, delimiter='\t'):
            if key in members or key not in intervals:
                raise ValueError('Invalid membership')
            members[key] = rep
    if set(members) != set(intervals):
        raise ValueError('Membership universe differs')
    remaining = dict(members)
    for key, rep in db.execute('SELECT interval_id,representative FROM members'):
        if remaining.pop(key,None) != rep or members.get(rep) != rep:
            raise ValueError('Stored membership differs')
    if remaining:
        raise ValueError('Missing membership')
    model_reps = defaultdict(set)
    models_by_rep = defaultdict(set)
    for key, (model, _, _) in intervals.items():
        model_reps[model].add(members[key])
        models_by_rep[members[key]].add(model)
    proteins = {}
    for gene, taxon, protein, model in source.execute('SELECT native_gene_id,taxon_id,protein_id,model FROM links'):
        if model in model_reps:
            proteins[gene] = (taxon,protein,model)
    remaining = dict(proteins)
    for gene,taxon,protein,model in db.execute('SELECT native_gene_id,taxon_id,protein_id,model FROM proteins'):
        if remaining.pop(gene,None) != (taxon,protein,model):
            raise ValueError('Protein source identity differs')
    if remaining or {v[2] for v in proteins.values()} != set(model_reps):
        raise ValueError('Missing protein/model links')
    cluster_genes = defaultdict(set)
    taxa = defaultdict(set)
    for gene,(taxon,_,model) in proteins.items():
        for rep in model_reps[model]:
            cluster_genes[rep].add(gene)
            taxa[rep].add(taxon)
    remaining = {(rep,gene) for rep,genes in cluster_genes.items() for gene in genes}
    expected_links = len(remaining)
    for row in db.execute('SELECT representative,native_gene_id FROM cluster_proteins'):
        if row not in remaining:
            raise ValueError('Unexpected or duplicate cluster/protein link')
        remaining.remove(row)
    if remaining:
        raise ValueError('Missing cluster/protein links')
    families = {}
    for guide in ('profile','mafft'):
        mapping = {gene:family for gene,family in source.execute('SELECT native_gene_id,family FROM family_links WHERE guide=?',(guide,)) if gene in proteins}
        if set(mapping) != set(proteins):
            raise ValueError('Family source incomplete')
        remaining = dict(mapping)
        for gene,family in db.execute('SELECT native_gene_id,family FROM family_links WHERE guide=?',(guide,)):
            if remaining.pop(gene,None) != family:
                raise ValueError('Family source mapping differs')
        if remaining:
            raise ValueError('Missing family mapping')
        families[guide] = mapping
    if db.execute('SELECT COUNT(*) FROM family_links').fetchone()[0] != 2*len(proteins):
        raise ValueError('Unexpected family guide')
    counts = Counter(members.values())
    remaining = {rep:(n,len(models_by_rep[rep])) for rep,n in counts.items()}
    for rep,n,models in db.execute('SELECT representative,intervals,models FROM cluster_intervals'):
        if remaining.pop(rep,None) != (n,models):
            raise ValueError('Stored interval summary differs')
    if remaining:
        raise ValueError('Missing interval summary')
    seen = set()
    with (root/'domain_cluster_composition.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rep = row['representative']
            if rep in seen or rep not in counts:
                raise ValueError('Unexpected summary cluster')
            genes = cluster_genes[rep]
            expected = [counts[rep],len(models_by_rep[rep]),len(genes),len(taxa[rep])]
            expected += [len({families[g][gene] for gene in genes}) for g in ('profile','mafft')]
            if [int(row[k]) for k in ('intervals','models','proteins','taxa','profile_families','mafft_families')] != expected:
                raise ValueError('Cluster composition differs')
            seen.add(rep)
    observed = dict(intervals=len(intervals),models=len(model_reps),proteins=len(proteins),clusters=len(seen),cluster_protein_links=expected_links)
    if seen != set(counts) or any(r[k] != v for k,v in observed.items()):
        raise ValueError('Receipt or summary universe differs')
    db.close()
    source.close()
    verify()
    verify_outputs()
    result = dict(status='passed_full_domain_cluster_composition_source_readback',plan_sha256=ph,producer_receipt_sha256=receipt_hash,**observed,
                  scope='Every interval, membership, selected source protein, both complete family mappings, cluster/protein links, database interval summary and TSV count reconstructed independently using Python sets and counters. Does not establish confidence qualification, homology, orthology or evolutionary events.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
