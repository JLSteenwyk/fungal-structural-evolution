#!/usr/bin/env python3
"""Independently compare all bridge identities and memberships to their sources."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def partitions(path):
    """Read the serialized one-family-per-line format without native MCL code."""
    with Path(path).open() as f:
        if f.readline().strip() != '(mclmatrix' or f.readline().strip() != 'begin':
            raise ValueError('Unexpected cluster header')
        index = 0
        for line in f:
            tokens = line.split()
            if tokens == [')']:
                if f.read().strip():
                    raise ValueError('Trailing cluster data')
                return
            if len(tokens) < 3 or tokens[0] != str(index) or tokens[-1] != '$':
                raise ValueError('Malformed cluster row')
            genes = tokens[1:-1]
            if len(genes) != len(set(genes)):
                raise ValueError('Duplicate cluster member')
            yield f'OG{index:07d}', sorted(genes)
            index += 1
        raise ValueError('Missing cluster terminator')


def audit(plan_path):
    plan = json.loads(plan_path.read_text())
    pins = {plan_path: sha(plan_path), **plan['pins']}
    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed audit input: '+str(path))
    verify()
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    start = time.time()
    def state(stage, **counts):
        p = out/'state.tmp'
        p.write_text(json.dumps(dict(stage=stage,elapsed_seconds=time.time()-start,**counts))+'\n')
        p.replace(out/'state.json')
    state('waiting_for_exact_producer')
    proc = Path('/proc')/str(plan['producer_pid'])/'stat'
    while proc.exists():
        try:
            fields = proc.read_text().rsplit(')',1)[1].split()
        except FileNotFoundError:
            break
        if fields[19] != str(plan['producer_start_ticks']) or fields[0] == 'Z':
            break
        time.sleep(20)
    source = json.loads(Path(plan['producer_plan']).read_text())
    receipt_path = Path(source['output'])/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    path = Path(source['output'])/'family_domain_bridge.sqlite'
    if (receipt['status'] != 'complete_family_domain_bridge_requires_independent_readback'
            or receipt['plan_sha256'] != sha(plan['producer_plan'])
            or receipt['bridge_sha256'] != sha(path)):
        raise ValueError('Incomplete or changed bridge')
    pins[receipt_path] = sha(receipt_path)
    pins[path] = receipt['bridge_sha256']
    for p, digest in source['pins'].items():
        pins[p] = digest
    verify()
    db = sqlite3.connect('file:'+str(path.resolve())+'?mode=ro',uri=True)
    raw = sqlite3.connect('file:'+str(Path(source['architecture_database']).resolve())+'?mode=ro',uri=True)
    species = {}
    for line in Path(source['species_ids']).read_text().splitlines():
        key, value = line.split(': ',1)
        species[key] = value.removesuffix('.faa')
    n = 0
    with open(source['sequence_ids']) as f:
        for line in f:
            gene, protein = line.rstrip('\n').split(': ',1)
            expected = (species[gene.split('_',1)[0]],protein)
            for table in ['native_proteins','proteins']:
                if db.execute(f'SELECT taxon_id,protein_id FROM {table} WHERE native_gene_id=?',(gene,)).fetchone() != expected:
                    raise ValueError('Changed native identity')
            n += 1
            if n % 100000 == 0:
                state('native_identity_readback', proteins=n)
    if n != source['proteins'] or len(set(species.values())) != source['taxa']:
        raise ValueError('Source scope mismatch')
    for table in ['native_proteins','proteins']:
        if db.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0] != n:
            raise ValueError('Protein scope mismatch')
    sql = 'SELECT taxon_id,protein_id,sequence_id,query_source FROM proteins ORDER BY taxon_id,protein_id'
    checked = 0
    for a,b in itertools.zip_longest(db.execute(sql),raw.execute(sql)):
        if a != b:
            raise ValueError('Changed architecture sequence/source link')
        checked += 1
    if checked != n:
        raise ValueError('Architecture scope mismatch')
    summaries = []
    for guide in source['guides']:
        families = proteins = 0
        with open(guide['crosswalk']) as f:
            for entry,row in itertools.zip_longest(partitions(guide['clusters']),csv.DictReader(f,delimiter='\t')):
                if entry is None or row is None:
                    raise ValueError('Family universe differs')
                name, genes = entry
                actual = [r[0] for r in db.execute('SELECT native_gene_id FROM assignments WHERE guide=? AND family=? ORDER BY native_gene_id',(guide['name'],name))]
                if actual != genes:
                    raise ValueError('Family membership differs')
                digest = hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest()
                expected = (guide['name'],name,row['source_type'],row['source_clade'],row['source_family'],len(genes),digest)
                if (row['new_family'] != name or int(row['proteins']) != len(genes)
                        or row['membership_sha256'] != digest
                        or db.execute('SELECT * FROM families WHERE guide=? AND family=?',(guide['name'],name)).fetchone() != expected):
                    raise ValueError('Family provenance differs')
                families += 1
                proteins += len(genes)
                if families % 10000 == 0:
                    state('family_readback',guide=guide['name'],families=families,proteins=proteins)
        if families != guide['families'] or proteins != n:
            raise ValueError('Partition scope differs')
        summaries.append(dict(guide=guide['name'],families=families,proteins=proteins))
    if (db.execute('SELECT COUNT(*) FROM assignments').fetchone()[0] != n*len(summaries)
            or db.execute('SELECT COUNT(*) FROM families').fetchone()[0] != sum(g['families'] for g in summaries)
            or db.execute('PRAGMA integrity_check').fetchone()[0] != 'ok'
            or db.execute('PRAGMA foreign_key_check').fetchone()):
        raise ValueError('Unexpected rows or database integrity failure')
    db.close(); raw.close()
    verify()
    result = dict(status='passed_complete_independent_family_domain_bridge_readback',proteins=n,taxa=source['taxa'],guides=summaries,
                  producer_receipt_sha256=sha(receipt_path),plan_sha256=sha(plan_path),script_sha256=sha(__file__),elapsed_seconds=time.time()-start,
                  scope='Every native identity, domain sequence/source link, family member and source crosswalk checked. Does not infer orthology or domain events.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    state(result['status'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    audit(p.parse_args().plan)
