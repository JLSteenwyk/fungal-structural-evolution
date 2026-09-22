#!/usr/bin/env python3
"""Independently validate every candidate architecture and every protein link."""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import time
import psutil
from audit_busco_gene_copies import sha

POLICIES = {'alignment_evalue', 'alignment_bitscore', 'envelope_evalue', 'envelope_bitscore'}
FIELDS = ['hit_id', 'pfam_accession', 'pfam_type', 'pfam_clan', 'alignment_start',
          'alignment_end', 'envelope_start', 'envelope_end', 'hmm_coverage']


def check_payload(payload, hits, policies):
    if set(payload) != {'alternatives', 'policies', 'status'} or set(payload['policies']) != POLICIES or set(policies) != POLICIES:
        raise ValueError('Payload or policy schema differs')
    expected_status = 'candidate_architecture_requires_validation' if hits else 'no_GA_hit_not_proven_absence'
    if payload['status'] != expected_status:
        raise ValueError('Missing-hit interpretation differs')
    by_id = {h['hit_id']: h for h in hits}
    if len(by_id) != len(hits):
        raise ValueError('Duplicate raw hit identity')
    signatures, references = set(), set()
    for index, alternative in enumerate(payload['alternatives']):
        if set(alternative) != {'annotations', 'ordered_model_tokens', 'model_multiplicities',
                                'token_signature_sha256', 'alignment_overlap_pairs',
                                'partial_hmm_hits_below_070', 'interpretation'}:
            raise ValueError('Alternative schema differs')
        annotated = alternative['annotations']
        ids = tuple(h['hit_id'] for h in annotated)
        if len(set(ids)) != len(ids) or not set(ids) <= set(by_id) or ids in signatures:
            raise ValueError('Duplicate, missing or repeated alternative identities')
        signatures.add(ids)
        expected = [{k: by_id[i][k] for k in FIELDS} for i in ids]
        for h in expected:
            for k in ['alignment_start', 'alignment_end', 'envelope_start', 'envelope_end']:
                h[k] = int(h[k])
        expected.sort(key=lambda h: (h['alignment_start'], h['alignment_end'], h['hit_id']))
        if annotated != expected:
            raise ValueError('Annotation fields or coordinate order differ')
        tokens = [[h['pfam_accession'], h['pfam_type']] for h in expected]
        multiplicities = Counter(tuple(t) for t in tokens)
        overlap = sum(max(a['alignment_start'], b['alignment_start']) <= min(a['alignment_end'], b['alignment_end'])
                      for a, b in itertools.combinations(expected, 2))
        partial = sum(Decimal(h['hmm_coverage']) < Decimal('0.70') for h in expected)
        token_hash = hashlib.sha256(json.dumps(tokens, separators=(',', ':'), sort_keys=True).encode()).hexdigest()
        if (alternative['ordered_model_tokens'] != tokens
                or alternative['model_multiplicities'] != [[a, t, n] for (a, t), n in sorted(multiplicities.items())]
                or alternative['alignment_overlap_pairs'] != overlap
                or alternative['partial_hmm_hits_below_070'] != partial
                or alternative['token_signature_sha256'] != token_hash
                or alternative['interpretation'] != 'Coordinate-sorted candidate annotations; repeated models and all Pfam types retained. Signature is not a homology assignment or validated biological architecture.'):
            raise ValueError('Candidate tokens, counts, overlaps or interpretation differ')
    for key in POLICIES:
        observed, source = payload['policies'][key], policies[key]
        index = observed['alternative_index']
        if not isinstance(index, int) or not 0 <= index < len(payload['alternatives']):
            raise ValueError('Invalid alternative reference')
        references.add(index)
        actual = [h['hit_id'] for h in payload['alternatives'][index]['annotations']]
        wanted = source['retained_hits']
        if len(wanted) != len(set(wanted)) or set(actual) != set(wanted):
            raise ValueError('Policy-specific retained membership differs')
        expected_state = dict(alternative_index=index,
                              primary_rank_tie_pairs=len(source['primary_rank_ties']),
                              unresolved_overlap_pairs=len(source['unresolved_overlap_pairs']),
                              candidate_nested_pairs=len(source['candidate_nested_pairs']))
        if observed != expected_state:
            raise ValueError('Policy uncertainty differs')
    if references != set(range(len(payload['alternatives']))):
        raise ValueError('Unused alternative')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    plan_sha = sha(a.plan)

    def verify():
        if sha(a.plan) != plan_sha:
            raise ValueError('Readback plan changed')
        for name, expected in plan['pins'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Pinned dependency changed: ' + name)

    verify()
    out = Path(plan['output'])
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    start = time.time()
    counts = Counter()

    def state(status):
        temp = out / 'state.tmp'
        temp.write_text(json.dumps(dict(status=status, counts=dict(counts), elapsed_seconds=time.time() - start)) + '\n')
        temp.replace(out / 'state.json')

    state('waiting_for_pinned_producer')
    while True:
        try:
            producer = psutil.Process(plan['producer_pid'])
            live = producer.create_time() == plan['producer_create_time'] and producer.status() != psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            live = False
        if not live:
            break
        time.sleep(20)
    verify()
    db_root, competition, architecture = [Path(plan[k]) for k in ['database', 'competition', 'architecture']]
    produced = json.loads((architecture / 'receipt.json').read_text())
    raw_receipt = json.loads((db_root / 'receipt.json').read_text())
    comp_receipt = json.loads((competition / 'receipt.json').read_text())
    comp_audit = json.loads((competition / 'readback.json').read_text())
    if (produced['status'] != 'completed_candidate_architecture_database_requires_independent_readback'
            or produced['database_receipt_sha256'] != sha(db_root / 'receipt.json')
            or produced['competition_readback_sha256'] != sha(competition / 'readback.json')
            or comp_audit['status'] != 'passed_complete_independent_domain_competition_readback'
            or comp_audit['production_receipt_sha256'] != sha(competition / 'receipt.json')):
        raise ValueError('Unbound production or competition receipts')
    required = {db_root / 'domains.sqlite': raw_receipt['artifacts']['domains.sqlite'],
                competition / 'query_competition.tsv.gz': comp_receipt['artifacts']['query_competition.tsv.gz'],
                architecture / 'candidate_architectures.sqlite': produced['artifacts']['candidate_architectures.sqlite'],
                architecture / 'receipt.json': sha(architecture / 'receipt.json')}
    for path, digest in required.items():
        if sha(path) != digest:
            raise ValueError('Changed data artifact')
    raw = sqlite3.connect('file:' + str((db_root / 'domains.sqlite').resolve()) + '?mode=ro', uri=True)
    derived = sqlite3.connect('file:' + str((architecture / 'candidate_architectures.sqlite').resolve()) + '?mode=ro', uri=True)
    raw.row_factory = sqlite3.Row
    sql = 'SELECT q.sequence_id,q.search_partition,' + ','.join('h.' + f for f in FIELDS) + ' FROM queries q LEFT JOIN hits h ON h.sequence_id=q.sequence_id ORDER BY q.sequence_id'
    groups = itertools.groupby(raw.execute(sql), key=lambda r: r['sequence_id'])
    observed = derived.execute('SELECT sequence_id,search_partition,raw_hits,retained_sets_agree,candidate_architectures_json FROM queries ORDER BY sequence_id')
    state('auditing_all_candidate_architectures')
    csv.field_size_limit(256 * 1024**2)
    with gzip.open(competition / 'query_competition.tsv.gz', 'rt') as handle:
        for group, record, saved in itertools.zip_longest(groups, csv.DictReader(handle, delimiter='\t'), observed):
            if group is None or record is None or saved is None:
                raise ValueError('Different query counts')
            sid, rows = group[0], list(group[1])
            hits = [dict(r) for r in rows if r['hit_id'] is not None]
            if (sid != record['sequence_id'] or sid != saved[0]
                    or rows[0]['search_partition'] != record['search_partition']
                    or saved[1] != record['search_partition']
                    or len(hits) != int(record['raw_hits']) or saved[2] != len(hits)
                    or saved[3] != int(record['retained_sets_agree'])):
                raise ValueError('Query identity, partition, hit count or agreement differs')
            payload = json.loads(saved[4])
            check_payload(payload, hits, json.loads(record['policies_json']))
            alternatives = len(payload['alternatives'])
            if (alternatives == 1) != bool(saved[3]):
                raise ValueError('Retained-set agreement differs')
            counts.update(queries=1, raw_hits=len(hits), alternatives=alternatives,
                          no_hit_queries=int(not hits), policy_disagreement_queries=int(not saved[3]))
            if counts['queries'] % 100000 == 0:
                state('auditing_all_candidate_architectures')
                print(dict(counts), flush=True)
    if dict(counts) != produced['counts'] or counts['queries'] != raw_receipt['queries']:
        raise ValueError('Production aggregate mismatch')
    state('auditing_all_protein_links')
    raw.row_factory = None
    query = 'SELECT taxon_id,protein_id,sequence_id,query_source FROM proteins ORDER BY taxon_id,protein_id'
    proteins, taxa = 0, set()
    for left, right in itertools.zip_longest(raw.execute(query), derived.execute(query)):
        if left is None or right is None or left != right:
            raise ValueError('Protein/taxon/sequence/source link differs')
        proteins += 1
        taxa.add(left[0])
    if proteins != produced['proteins'] or proteins != raw_receipt['proteins'] or len(taxa) != produced['taxa']:
        raise ValueError('Protein scope mismatch')
    if derived.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or derived.execute('PRAGMA foreign_key_check').fetchone():
        raise ValueError('Database integrity failure')
    raw.close()
    derived.close()
    verify()
    for path, digest in required.items():
        if sha(path) != digest:
            raise ValueError('Artifact changed during readback')
    result = dict(status='passed_complete_independent_candidate_architecture_readback',
                  counts=dict(counts), proteins=proteins, taxa=len(taxa),
                  production_receipt_sha256=sha(architecture / 'receipt.json'),
                  plan_sha256=plan_sha, script_sha256=sha(Path(__file__)),
                  elapsed_seconds=time.time() - start,
                  scope='Every query alternative, annotation field, coordinate order, repeat/type count, signature, uncertainty flag and policy assignment independently checked; every protein/taxon link matched. Biological architecture and gains/losses remain unvalidated.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    state(result['status'])


if __name__ == '__main__':
    main()
