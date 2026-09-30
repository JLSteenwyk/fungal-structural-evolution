#!/usr/bin/env python3
"""Independently rebuild full context/measurement work states using SQL gene-identity joins."""
import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from reference_context_measurement_sources import load_sources
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def digest(ends):
    return hashlib.sha256(json.dumps(sorted(ends), separators=(',', ':')).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    source, queue, inventory, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    out = Path(plan['output']); r = json.loads((out / 'receipt.json').read_text())
    assert r['status'] == 'complete_full_reference_context_measurement_design_pending_independent_readback' and r['plan_sha256'] == sha(args.plan)
    assert r['source_hashes'] == bindings and set(r['artifacts']) == {'context_measurement_design.jsonl'}
    bind(bindings, out / 'receipt.json')
    for name, h in r['artifacts'].items(): bind(bindings, out / name, h)
    verify(bindings)
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE duplicates (guide TEXT,family TEXT,node TEXT,lo TEXT,hi TEXT,payload TEXT,PRIMARY KEY(guide,family,node,lo,hi))')
    db.execute('CREATE TABLE pairs (key TEXT PRIMARY KEY,a TEXT,av INT,b TEXT,bv INT,work TEXT)')
    db.execute('CREATE TABLE links (guide TEXT,family TEXT,node TEXT,lo TEXT,hi TEXT,ref TEXT,focal TEXT,payload TEXT,checked INT DEFAULT 0,PRIMARY KEY(guide,family,node,lo,hi,ref,focal))')
    db.execute('CREATE TABLE seen (guide TEXT,ordinal INT,family TEXT,node TEXT,lo TEXT,hi TEXT,PRIMARY KEY(guide,ordinal),UNIQUE(guide,family,node,lo,hi))')
    with (queue / 'event_model_pair_links.tsv').open() as handle:
        db.executemany('INSERT INTO duplicates VALUES (?,?,?,?,?,?)',
                       ((x['guide'], x['family'], x['gene_node'], *sorted([x['gene_a'], x['gene_b']]), json.dumps(x)) for x in csv.DictReader(handle, delimiter='\t')))
    with (inventory / 'model_pairs.tsv').open() as handle:
        for x in csv.DictReader(handle, delimiter='\t'):
            ends = [(x['model_a'], int(x['version_a'])), (x['model_b'], int(x['version_b']))]
            assert ends[0] != ends[1] and digest(ends) == x['pair_key']
            db.execute('INSERT INTO pairs VALUES (?,?,?,?,?,?)', (x['pair_key'], *ends[0], *ends[1], x['work_disposition']))
    with (inventory / 'event_reference_comparisons.tsv').open() as handle:
        db.executemany('INSERT INTO links(guide,family,node,lo,hi,ref,focal,payload) VALUES (?,?,?,?,?,?,?,?)',
                       ((x['guide'], x['family'], x['gene_node'], *sorted([x['gene_a'], x['gene_b']]), x['reference_gene'], x['focal_gene'], json.dumps(x)) for x in csv.DictReader(handle, delimiter='\t')))
    assert db.execute('SELECT COUNT(*) FROM duplicates').fetchone()[0] == plan['expected']['target_contexts']
    assert db.execute('SELECT COUNT(*) FROM pairs').fetchone()[0] == plan['expected']['model_pairs']
    assert db.execute('SELECT COUNT(*) FROM links').fetchone()[0] == plan['expected']['availability_side_links']
    guide_counts, counts = Counter(), Counter(); contexts = ties = logical = 0
    with source.open() as original, (out / 'context_measurement_design.jsonl').open() as output:
        for left, right in zip_longest(original, output):
            assert left is not None and right is not None
            context = json.loads(left); actual = json.loads(right); row = context['source']; guide = context['source_guide']
            assert set(actual) == {'native_context', 'duplicate_models', 'duplicate_pair_key', 'duplicate_comparison_status', 'measurement_designs'}
            assert actual['native_context'] == context and context['source_row_sha256'] == hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            prefix = (guide, row['family'], row['gene_node'], *sorted([row['gene_a'], row['gene_b']]))
            q = db.execute('SELECT payload FROM duplicates WHERE guide=? AND family=? AND node=? AND lo=? AND hi=?', prefix).fetchone(); assert q
            queued = json.loads(q[0]); mapped = {}
            for side in ['a', 'b']:
                matches = [s for s in ['a', 'b'] if queued['gene_' + s] == row['gene_' + side]]; assert len(matches) == 1
                qside = matches[0]
                assert queued['model_' + qside] == row['model_' + side]
                mapped[side] = dict(gene=row['gene_' + side], model_id=queued['model_' + qside], version=int(queued['version_' + qside]))
            assert queued['taxon_id'] == row['taxon_id']
            dupends = [(mapped[s]['model_id'], mapped[s]['version']) for s in ['a', 'b']]
            assert digest(dupends) == queued['pair_key']
            assert bool(int(queued['same_model'])) == (dupends[0] == dupends[1])
            assert actual['duplicate_models'] == mapped and actual['duplicate_pair_key'] == queued['pair_key'] and actual['duplicate_comparison_status'] == queued['comparison_status']
            eligible = row['status'] in ['provisional_reference_available', 'no_modeled_nonfocal_sister']
            assert context['parent_context_eligible'] == eligible
            assert set(actual['measurement_designs']) == {'availability', 'sequence_first'}
            for design in ['availability', 'sequence_first']:
                refs = context['designs'][design]; exported = actual['measurement_designs'][design]
                assert len(refs) == len(exported) and [x['reference_gene'] for x in refs] == sorted(set(x['reference_gene'] for x in refs))
                if not refs: counts[f'{guide}|{design}|parent={int(eligible)}|no_reference_gene'] += 1
                for index, (ref, work) in enumerate(zip(refs, exported)):
                    assert ref['lexical_choice'] == (index == 0) and set(work) == {'reference', 'side_work'} and work['reference'] == ref
                    assert set(work['side_work']) == {'a', 'b'}
                    assert all(v in ['both', 'only_a', 'only_b', 'neither'] for v in ref['native_coorthology'].values())
                    for side in ['a', 'b']:
                        focal = mapped[side]['model_id'], mapped[side]['version']
                        reference = (ref['reference_model'], int(ref['reference_version'])) if ref['reference_model'] else None
                        if not reference: assert ref['reference_version'] == ''
                        key = digest([focal, reference]) if reference else ''
                        sql = db.execute('SELECT a,av,b,bv,work FROM pairs WHERE key=?', (key,)).fetchone() if key else None
                        endpoint, ledger_work = '', ''
                        if reference is None: disposition = 'reference_model_missing'
                        elif reference == focal: disposition = 'identical_model_not_independent'
                        elif sql:
                            assert sorted([focal, reference]) == sorted([(sql[0], sql[1]), (sql[2], sql[3])])
                            endpoint = 0 if focal == (sql[0], sql[1]) else 1; ledger_work = sql[4]
                            disposition = 'in_full_reference_measurement_design'
                        else: disposition = 'outside_full_reference_measurement_design'
                        lk = (*prefix, ref['reference_gene'], mapped[side]['gene'])
                        found = db.execute('SELECT payload,checked FROM links WHERE guide=? AND family=? AND node=? AND lo=? AND hi=? AND ref=? AND focal=?', lk).fetchone() if design == 'availability' else None
                        if found:
                            link = json.loads(found[0]); assert found[1] == 0
                            assert (link['focal_model'], int(link['focal_version'])) == focal and reference is not None and (link['reference_model'], int(link['reference_version'])) == reference
                            assert link['focal_gene'] == link['gene_' + link['focal_side']] and link['pair_key'] == key
                            assert int(link['lexical_representative']) == int(ref['lexical_choice'])
                            assert link['work_disposition'] == ('identical_model' if focal == reference else ledger_work)
                            db.execute('UPDATE links SET checked=1 WHERE guide=? AND family=? AND node=? AND lo=? AND hi=? AND ref=? AND focal=?', lk)
                        if design == 'availability' and eligible: assert found and disposition in ['in_full_reference_measurement_design', 'identical_model_not_independent']
                        expected = dict(duplicate_gene=mapped[side]['gene'], duplicate_model=focal[0], duplicate_version=focal[1], reference_model=reference[0] if reference else '',
                                        reference_version=reference[1] if reference else '', pair_key=key, measurement_disposition=disposition, current_pair_focal_endpoint=endpoint,
                                        pair_work_disposition=ledger_work, availability_ledger_link=bool(found))
                        assert work['side_work'][side] == expected
                        counts[f'{guide}|{design}|parent={int(eligible)}|{disposition}'] += 1; logical += 1
                    ties += 1
            db.execute('INSERT INTO seen VALUES (?,?,?,?,?,?)', (guide, context['source_row_number'], *prefix[1:]))
            guide_counts[guide] += 1; contexts += 1
            if contexts % 25000 == 0: print('Independent SQL reference context/model design', contexts, '/', plan['expected']['target_contexts'], flush=True)
    assert db.execute('SELECT COUNT(*) FROM links WHERE checked=1').fetchone()[0] == plan['expected']['availability_side_links']
    assert db.execute('SELECT COUNT(*) FROM duplicates d LEFT JOIN seen s USING(guide,family,node,lo,hi) WHERE s.ordinal IS NULL').fetchone()[0] == 0
    assert dict(guide_counts) == upstream['guide_contexts'] and contexts == plan['expected']['target_contexts']
    for guide, n in guide_counts.items(): assert db.execute('SELECT MIN(ordinal),MAX(ordinal),COUNT(*) FROM seen WHERE guide=?', (guide,)).fetchone() == (1, n, n)
    assert ties == plan['expected']['reference_tie_records'] and logical == plan['expected']['duplicate_reference_links']
    summary = dict(target_contexts=contexts, context_design_records=2 * contexts, reference_tie_records=ties, duplicate_reference_links=logical,
                   availability_side_links_checked=plan['expected']['availability_side_links'], guide_contexts=dict(guide_counts), measurement_disposition_counts=dict(counts))
    assert all(r[k] == v for k, v in summary.items()); verify(bindings)
    result = dict(status='passed_full_reference_context_measurement_design_sql_readback', plan_sha256=sha(args.plan), producer_receipt_sha256=sha(out / 'receipt.json'),
                  checker_sha256=sha(__file__), **summary, source_hashes=bindings, scientific_eligibility=False,
                  scope='All source contexts/native fields preserved and every duplicate gene/model/version joined independently with uniquely keyed SQL queue records. Every design/tie/lexical/parent/missing-model state, pair hash/current endpoint, full-ledger disposition and all availability gene-side links independently reconstructed. Both full context and availability-ledger universes exhausted, all row ordinals and aggregates checked. Shares source/proof I/O only; no producer projection functions. No measurement, biological orthology or asymmetry acceptance.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
