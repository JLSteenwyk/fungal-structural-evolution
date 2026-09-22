#!/usr/bin/env python3
"""Preserve native family membership and taxon-specific candidate-domain links."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sqlite3
import time
from orthofinder.tools import mcl


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def build(plan_path):
    plan = json.loads(plan_path.read_text())
    plan_hash = sha(plan_path)

    def verify():
        if sha(plan_path) != plan_hash:
            raise ValueError('Plan changed')
        for path, digest in plan['pins'].items():
            if sha(path) != digest:
                raise ValueError('Changed input: ' + path)

    verify()
    merge_audit = json.loads(Path(plan['merge_readback']).read_text())
    if (merge_audit['status'] != 'passed_complete_guide_discovery_merge_readback'
            or merge_audit['merge_receipt_sha256'] != sha(plan['merge_receipt'])):
        raise ValueError('Unverified family partition')
    audit = json.loads(Path(plan['architecture_readback']).read_text())
    receipt = json.loads(Path(plan['architecture_receipt']).read_text())
    if (audit['status'] != 'passed_complete_independent_candidate_architecture_readback'
            or audit['production_receipt_sha256'] != sha(plan['architecture_receipt'])
            or sha(plan['architecture_database']) != receipt['artifacts']['candidate_architectures.sqlite']):
        raise ValueError('Unverified candidate architectures')
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient disk')
    out.mkdir()
    start = time.time()

    def state(stage, **counts):
        p = out / 'state.tmp'
        p.write_text(json.dumps(dict(stage=stage, elapsed_seconds=time.time()-start, **counts))+'\n')
        p.replace(out / 'state.json')

    db = sqlite3.connect(out / 'family_domain_bridge.sqlite', uri=True)
    db.execute('PRAGMA foreign_keys=ON')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('PRAGMA cache_size=-65536')
    db.execute('ATTACH DATABASE ? AS architecture',
               ('file:' + str(Path(plan['architecture_database']).resolve()) + '?mode=ro',))
    db.execute('CREATE TABLE native_proteins(native_gene_id TEXT PRIMARY KEY,taxon_id TEXT NOT NULL,protein_id TEXT NOT NULL,UNIQUE(taxon_id,protein_id))')
    species = {}
    for line in Path(plan['species_ids']).read_text().splitlines():
        key, value = line.split(': ', 1)
        if key in species or not value.endswith('.faa'):
            raise ValueError('Malformed/duplicate species identifier')
        species[key] = value[:-4]
    if len(species) != plan['taxa'] or len(set(species.values())) != plan['taxa']:
        raise ValueError('Species scope differs')
    state('native_protein_mapping')
    n = 0
    with open(plan['sequence_ids']) as f:
        for line in f:
            gene, protein = line.rstrip('\n').split(': ', 1)
            db.execute('INSERT INTO native_proteins VALUES (?,?,?)',
                       (gene, species[gene.split('_', 1)[0]], protein))
            n += 1
            if n % 100000 == 0:
                db.commit()
                state('native_protein_mapping', proteins=n)
    db.commit()
    if n != plan['proteins']:
        raise ValueError('Native protein scope differs')
    db.execute('CREATE TABLE proteins(native_gene_id TEXT PRIMARY KEY REFERENCES native_proteins(native_gene_id),taxon_id TEXT NOT NULL,protein_id TEXT NOT NULL,sequence_id TEXT NOT NULL,query_source TEXT NOT NULL,UNIQUE(taxon_id,protein_id))')
    db.execute('INSERT INTO proteins SELECT n.native_gene_id,n.taxon_id,n.protein_id,a.sequence_id,a.query_source FROM native_proteins n JOIN architecture.proteins a USING(taxon_id,protein_id)')
    if (db.execute('SELECT COUNT(*) FROM proteins').fetchone()[0] != n
            or db.execute('SELECT COUNT(*) FROM architecture.proteins').fetchone()[0] != n):
        raise ValueError('Native and domain protein universes differ')
    db.execute('CREATE TABLE families(guide TEXT,family TEXT,source_type TEXT NOT NULL,source_clade TEXT NOT NULL,source_family TEXT NOT NULL,proteins INTEGER NOT NULL,membership_sha256 TEXT NOT NULL,PRIMARY KEY(guide,family))')
    db.execute('CREATE TABLE assignments(guide TEXT,native_gene_id TEXT REFERENCES proteins(native_gene_id),family TEXT NOT NULL,PRIMARY KEY(guide,native_gene_id),FOREIGN KEY(guide,family) REFERENCES families(guide,family))')
    summaries = []
    for guide in plan['guides']:
        state('family_assignments', guide=guide['name'])
        groups = mcl.GetPredictedOGs(guide['clusters'])
        total = 0
        with open(guide['crosswalk']) as f:
            rows = csv.DictReader(f, delimiter='\t')
            for index, pair in enumerate(itertools.zip_longest(groups, rows)):
                genes, row = pair
                if genes is None or row is None:
                    raise ValueError('Crosswalk and partition sizes differ')
                family = f'OG{index:07d}'
                digest = hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest()
                if row['new_family'] != family or int(row['proteins']) != len(genes) or row['membership_sha256'] != digest:
                    raise ValueError('Family crosswalk mismatch')
                db.execute('INSERT INTO families VALUES (?,?,?,?,?,?,?)',
                           (guide['name'], family, row['source_type'], row['source_clade'], row['source_family'], len(genes), digest))
                db.executemany('INSERT INTO assignments VALUES (?,?,?)',
                               ((guide['name'], gene, family) for gene in genes))
                total += len(genes)
                if index % 10000 == 0:
                    db.commit()
                    state('family_assignments', guide=guide['name'], families=index+1, proteins=total)
        if len(groups) != guide['families'] or total != n:
            raise ValueError('Partition scope differs')
        summaries.append(dict(guide=guide['name'], families=len(groups), proteins=total))
        del groups
        db.commit()
    state('indexing_and_integrity')
    db.execute('CREATE INDEX assignments_family ON assignments(guide,family)')
    db.execute('CREATE INDEX proteins_sequence ON proteins(sequence_id)')
    db.commit()
    if db.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or db.execute('PRAGMA foreign_key_check').fetchone():
        raise ValueError('Integrity failure')
    db.close()
    verify()
    result = dict(status='complete_family_domain_bridge_requires_independent_readback',
                  proteins=n, taxa=len(species), guides=summaries, plan_sha256=plan_hash,
                  script_sha256=sha(__file__), architecture_database_sha256=sha(plan['architecture_database']),
                  bridge_sha256=sha(out/'family_domain_bridge.sqlite'), elapsed_seconds=time.time()-start,
                  scope='Exact taxon/protein/sequence links and both complete family partitions. Attach the pinned architecture database to retrieve alternatives. Family assignments are not reconciled orthology; domain events remain uninferred.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    state(result['status'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    build(p.parse_args().plan)
