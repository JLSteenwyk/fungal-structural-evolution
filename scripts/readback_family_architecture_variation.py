#!/usr/bin/env python3
"""Reconstruct every family/policy summary directly from audited annotations."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sqlite3
import time

POLICIES = ('alignment_evalue', 'alignment_bitscore', 'envelope_evalue', 'envelope_bitscore')
FLAGS = ('rank_tie_proteins', 'unresolved_overlap_proteins', 'candidate_nested_proteins',
         'alignment_overlap_proteins', 'partial_hmm_proteins')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def variation(records, prefix):
    # Sorted token tuples preserve multiplicity without using the producer's
    # token-count representation. Records contain only annotated proteins.
    orders = {tokens for _, tokens in records}
    bags = {tuple(sorted(tokens)) for tokens in orders}
    kinds = {frozenset(tokens) for tokens in orders}
    by_bag = Counter(tuple(sorted(tokens)) for tokens in orders)
    by_kind = Counter(frozenset(bag) for bag in bags)
    taxon_signatures = {}
    for taxon, tokens in records:
        taxon_signatures.setdefault(taxon, set()).add(tokens)
    return {prefix+k: v for k, v in dict(
        proteins=len(records), taxa=len(taxon_signatures), ordered_signatures=len(orders),
        multisets=len(bags), model_type_sets=len(kinds),
        multisets_with_order_variation=sum(n > 1 for n in by_bag.values()),
        model_type_sets_with_multiplicity_variation=sum(n > 1 for n in by_kind.values()),
        taxa_with_multiple_signatures=sum(len(s) > 1 for s in taxon_signatures.values())).items()}


def reconstruct(source_rows):
    records = [(taxon, seq, int(hits), int(agree), json.loads(payload))
               for taxon, seq, hits, agree, payload in source_rows]
    copies = Counter(r[0] for r in records)
    result = []
    for policy in POLICIES:
        flags = Counter()
        observed, conservative = [], []
        for taxon, seq, hits, agree, payload in records:
            if not hits:
                continue
            entry = payload['policies'][policy]
            alt = payload['alternatives'][entry['alternative_index']]
            tokens = tuple(map(tuple, alt['ordered_model_tokens']))
            if not tokens:
                raise ValueError('Annotated query has no retained tokens')
            observed.append((taxon, tokens))
            values = [entry['primary_rank_tie_pairs'], entry['unresolved_overlap_pairs'],
                      entry['candidate_nested_pairs'], alt['alignment_overlap_pairs'],
                      alt['partial_hmm_hits_below_070']]
            for key, value in zip(FLAGS, values):
                flags[key] += value > 0
            if agree and not any(value > 0 for value in values):
                conservative.append((taxon, tokens))
        row = dict(policy=policy, proteins=len(records), taxa=len(copies),
                   unique_sequences=len({r[1] for r in records}),
                   taxa_with_multiple_family_members=sum(n > 1 for n in copies.values()),
                   maximum_members_per_taxon=max(copies.values(), default=0),
                   no_ga_hit_proteins=sum(r[2] == 0 for r in records),
                   policy_disagreement_proteins=sum(not r[3] for r in records),
                   **{key: flags[key] for key in FLAGS})
        row.update(variation(observed, 'observed_'))
        row.update(variation(conservative, 'conservative_'))
        result.append(row)
    return result


def check_family(table, guide, family, source_rows):
    metrics = reconstruct(source_rows)
    for expected in metrics:
        expected = dict(guide=guide, family=family, **expected)
        actual = next(table, None)
        if actual != {k: str(v) for k, v in expected.items()}:
            raise ValueError('Family/policy metrics differ: ' + guide + '/' + family)
    return metrics[0]['proteins']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', required=True, type=Path)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    pins = {str(a.plan): sha(a.plan), **plan['pins']}
    def verify():
        for path, digest in pins.items():
            if sha(path) != digest:
                raise ValueError('Changed source: ' + path)
    verify()
    out = Path(plan['output'])
    out.mkdir(exist_ok=False)
    started = time.time()
    def state(stage, **details):
        p = out/'state.tmp'
        p.write_text(json.dumps(dict(stage=stage, elapsed_seconds=time.time()-started, **details))+'\n')
        p.replace(out/'state.json')
    state('waiting_for_exact_variation_producer')
    proc = Path('/proc')/str(plan['producer_pid'])/'stat'
    while proc.exists():
        try:
            fields = proc.read_text().rsplit(')', 1)[1].split()
        except FileNotFoundError:
            break
        if fields[19] != str(plan['producer_start_ticks']) or fields[0] == 'Z':
            break
        time.sleep(20)
    source = json.loads(Path(plan['producer_plan']).read_text())
    folder = Path(source['output'])
    receipt_path = folder/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    if receipt['status'] != 'complete_family_architecture_variation_inventory':
        raise ValueError('Incomplete producer')
    if receipt['input_hashes'].get(plan['producer_plan']) != sha(plan['producer_plan']):
        raise ValueError('Producer plan mismatch')
    pins.update(receipt['input_hashes'])
    pins[str(receipt_path)] = sha(receipt_path)
    pins.update({str(folder/name): digest for name, digest in receipt['artifacts'].items()})
    verify()
    if shutil.disk_usage(out).free < plan['resources']['minimum_free_disk_gib']*2**30:
        raise ValueError('Insufficient free disk')
    bridge = json.loads(Path(source['bridge_receipt']).read_text())
    bridge_audit = json.loads(Path(source['bridge_readback']).read_text())
    if (bridge_audit['status'] != 'passed_complete_independent_family_domain_bridge_readback'
            or bridge_audit['producer_receipt_sha256'] != sha(source['bridge_receipt'])
            or bridge['bridge_sha256'] != sha(source['bridge_database'])
            or bridge['architecture_database_sha256'] != sha(source['architecture_database'])):
        raise ValueError('Unverified bridge inputs')
    db = sqlite3.connect('file:'+str(Path(source['bridge_database']).resolve())+'?mode=ro', uri=True)
    db.execute('ATTACH DATABASE ? AS architecture',
               ('file:'+str(Path(source['architecture_database']).resolve())+'?mode=ro',))
    summaries = []
    for guide in bridge['guides']:
        name = guide['guide']
        state('full_family_readback', guide=name, families=0)
        family_count = protein_count = 0
        # The bridge's membership itself has already passed a full independent
        # source readback; this audit checks the annotation summary on that frame.
        sql = '''SELECT a.family,p.taxon_id,p.sequence_id,q.raw_hits,
                        q.retained_sets_agree,q.candidate_architectures_json
                 FROM assignments a INDEXED BY assignments_family
                 JOIN proteins p ON a.native_gene_id=p.native_gene_id
                 JOIN architecture.queries q ON p.sequence_id=q.sequence_id
                 WHERE a.guide=? ORDER BY a.family'''
        with gzip.open(folder/(name+'_family_architectures.tsv.gz'), 'rt') as f:
            table = iter(csv.DictReader(f, delimiter='\t'))
            for family, group in itertools.groupby(db.execute(sql, (name,)), lambda r: r[0]):
                protein_count += check_family(table, name, family, (r[1:] for r in group))
                family_count += 1
                if family_count % 10000 == 0:
                    state('full_family_readback', guide=name, families=family_count, proteins=protein_count)
            if next(table, None) is not None:
                raise ValueError('Extra family/policy rows')
        summary = dict(guide=name, families=family_count, proteins=protein_count, rows=family_count*4)
        expected = next(r for r in receipt['guides'] if r['guide'] == name)
        if (summary != expected or family_count != guide['families']
                or protein_count != bridge['proteins']):
            raise ValueError('Incomplete family universe')
        summaries.append(summary)
    db.close()
    verify()
    result = dict(status='passed_full_family_architecture_variation_readback', guides=summaries,
                  producer_receipt_sha256=sha(receipt_path), plan_sha256=sha(a.plan),
                  script_sha256=sha(__file__), elapsed_seconds=time.time()-started,
                  scope='Every family/policy row reconstructed from the audited bridge and candidate annotation database without importing producer metrics. Uses the same SQLite input join; does not repeat Pfam search, membership inference or architecture resolution. No evolutionary event or significance claim.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    state(result['status'])


if __name__ == '__main__':
    main()
