#!/usr/bin/env python3
"""Join independently checked competition to raw annotations and every protein link."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import shutil
import sqlite3
import time
from candidate_domain_architecture import compact, describe
from audit_busco_gene_copies import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    plan_sha = sha(a.plan)

    def verify():
        if sha(a.plan) != plan_sha:
            raise ValueError('Plan changed')
        for name, expected in plan['pins'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Changed input: ' + name)

    verify()
    database, competition, out = [Path(plan[k]) for k in ['database', 'competition', 'output']]
    receipt = json.loads((competition / 'receipt.json').read_text())
    audit = json.loads((competition / 'readback.json').read_text())
    db_receipt = json.loads((database / 'receipt.json').read_text())
    db_audit = json.loads((database / 'readback.json').read_text())
    if (db_audit['status'] != 'passed_complete_domain_database_source_readback'
            or db_audit['database_receipt_sha256'] != sha(database / 'receipt.json')):
        raise ValueError('Domain database lacks bound source readback')
    if (audit['status'] != 'passed_complete_independent_domain_competition_readback'
            or audit['production_receipt_sha256'] != sha(competition / 'receipt.json')
            or audit['totals'] != receipt['totals']):
        raise ValueError('Competition lacks bound independent readback')
    if sha(database / 'domains.sqlite') != db_receipt['artifacts']['domains.sqlite']:
        raise ValueError('Domain database changed')
    if sha(competition / 'query_competition.tsv.gz') != receipt['artifacts']['query_competition.tsv.gz']:
        raise ValueError('Competition output changed')
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient disk')
    out.mkdir()
    start = time.time()
    source = sqlite3.connect('file:' + str((database / 'domains.sqlite').resolve()) + '?mode=ro', uri=True)
    source.row_factory = sqlite3.Row
    target = out / 'candidate_architectures.sqlite'
    dest = sqlite3.connect(target)
    dest.execute('PRAGMA foreign_keys=ON')
    dest.execute('CREATE TABLE queries(sequence_id TEXT PRIMARY KEY, search_partition TEXT NOT NULL, raw_hits INTEGER NOT NULL, retained_sets_agree INTEGER NOT NULL, candidate_architectures_json TEXT NOT NULL)')
    dest.execute('CREATE TABLE proteins(taxon_id TEXT NOT NULL, protein_id TEXT NOT NULL, sequence_id TEXT NOT NULL REFERENCES queries(sequence_id), query_source TEXT NOT NULL, PRIMARY KEY(taxon_id,protein_id))')
    fields = ['hit_id', 'pfam_accession', 'pfam_type', 'pfam_clan', 'alignment_start',
              'alignment_end', 'envelope_start', 'envelope_end', 'hmm_coverage']
    cursor = source.execute('SELECT q.sequence_id,q.search_partition,' + ','.join('h.' + k for k in fields) +
                            ' FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id')
    groups = itertools.groupby(cursor, key=lambda r: r['sequence_id'])
    counts = Counter()
    csv.field_size_limit(256 * 1024**2)

    def state(status):
        temp = out / 'state.tmp'
        temp.write_text(json.dumps(dict(status=status, counts=dict(counts), elapsed_seconds=time.time() - start)) + '\n')
        temp.replace(out / 'state.json')

    with gzip.open(competition / 'query_competition.tsv.gz', 'rt') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        for group, row in itertools.zip_longest(groups, rows):
            if group is None or row is None or group[0] != row['sequence_id']:
                raise ValueError('Query universes or order differ')
            sid, data = group[0], list(group[1])
            hits = [dict(r) for r in data if r['hit_id'] is not None]
            if len(hits) != int(row['raw_hits']) or data[0]['search_partition'] != row['search_partition']:
                raise ValueError('Hit count or partition differs')
            for hit in hits:
                for field in ['alignment_start', 'alignment_end', 'envelope_start', 'envelope_end']:
                    hit[field] = int(hit[field])
            payload = describe(hits, json.loads(row['policies_json']))
            agree = len(payload['alternatives']) == 1
            if agree != bool(int(row['retained_sets_agree'])):
                raise ValueError('Policy agreement differs')
            dest.execute('INSERT INTO queries VALUES (?,?,?,?,?)',
                         (sid, row['search_partition'], len(hits), int(agree), compact(payload)))
            counts.update(queries=1, raw_hits=len(hits), alternatives=len(payload['alternatives']),
                          no_hit_queries=int(not hits), policy_disagreement_queries=int(not agree))
            if counts['queries'] % 100000 == 0:
                dest.commit()
                state('building_candidate_architectures')
                print(dict(counts), flush=True)
                if target.stat().st_size > plan['resources']['output_allowance_gib'] * 2**30:
                    raise ValueError('Output allowance exceeded')
                if shutil.disk_usage(out).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
                    raise ValueError('Disk gate reached')
    dest.commit()
    if counts['queries'] != db_receipt['queries'] or counts['raw_hits'] != db_receipt['total_hit_rows']:
        raise ValueError('Incomplete annotation scope')
    state('copying_all_protein_links')
    source.row_factory = None
    dest.executemany('INSERT INTO proteins VALUES (?,?,?,?)',
                     source.execute('SELECT taxon_id,protein_id,sequence_id,query_source FROM proteins'))
    dest.execute('CREATE INDEX proteins_sequence ON proteins(sequence_id)')
    dest.commit()
    nproteins, ntaxa = dest.execute('SELECT COUNT(*),COUNT(DISTINCT taxon_id) FROM proteins').fetchone()
    if nproteins != db_receipt['proteins'] or ntaxa != db_receipt['taxa']:
        raise ValueError('Incomplete protein scope')
    if dest.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or dest.execute('PRAGMA foreign_key_check').fetchone():
        raise ValueError('Database integrity failure')
    dest.close()
    source.close()
    verify()
    result = dict(status='completed_candidate_architecture_database_requires_independent_readback',
                  counts=dict(counts), proteins=nproteins, taxa=ntaxa,
                  database_receipt_sha256=sha(database / 'receipt.json'),
                  competition_readback_sha256=sha(competition / 'readback.json'),
                  plan_sha256=plan_sha, script_sha256=sha(Path(__file__)),
                  artifacts={target.name: sha(target)}, elapsed_seconds=time.time() - start,
                  scope='All four policy-specific candidate architectures and all protein links. Raw model types, repeats, coordinates and uncertainty retained. No-hit is not absence; ordering and signatures do not establish validated architectures, homology, gains or losses.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    state(result['status'])


if __name__ == '__main__':
    main()
