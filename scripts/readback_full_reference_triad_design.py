#!/usr/bin/env python3
"""Independently rebuild every full triad and link using keyed SQL pair catalogs."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from reference_triad_design_sources import load_sources, SUMMARY_FIELDS
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def hash_value(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def sql_edge(db, left, right, roles, kind):
    key = hash_value(sorted([left, right])); row = db.execute('SELECT a,av,b,bv,work FROM pairs WHERE kind=? AND key=?', (kind, key)).fetchone()
    if left == right:
        status = 'identical_model_not_independent'; endpoints = direction = work = None
    elif row is None:
        status = 'outside_full_' + kind + '_measurement_design'; endpoints = direction = work = None
    else:
        endpoints = [[row[0], row[1]], [row[2], row[3]]]
        assert sorted(endpoints) == sorted([list(left), list(right)])
        direction = int(endpoints[0] != list(left)); work = row[4]
        status = 'in_full_' + kind + '_measurement_design'
    return dict(pair_key=key, measurement_disposition=status, desired_role_order=roles,
                current_pair_endpoints=endpoints, desired_direction_order=direction, pair_work_disposition=work)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    source, queue, inventory, upstream, bindings = load_sources(plan); bind(bindings, args.plan)
    out = Path(plan['output']); rp = out / 'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_reference_triad_work_design_pending_independent_readback' and receipt['plan_sha256'] == sha(args.plan)
    assert receipt['source_hashes'] == bindings and set(receipt['artifacts']) == {'context_triad_design.jsonl.gz', 'ordered_model_triads.jsonl'}
    bind(bindings, rp)
    for name, value in receipt['artifacts'].items(): bind(bindings, out / name, value)
    verify(bindings)
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE pairs(kind TEXT,key TEXT,a TEXT,av INT,b TEXT,bv INT,work TEXT,PRIMARY KEY(kind,key))')
    for kind, root in [('primary', queue), ('reference', inventory)]:
        with (root / 'model_pairs.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                ends = [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]
                assert ends[0] != ends[1] and row['pair_key'] == hash_value(sorted(ends))
                db.execute('INSERT INTO pairs VALUES (?,?,?,?,?,?,?)', (kind, row['pair_key'], *ends[0], *ends[1], row.get('work_disposition', 'primary_pair')))
        assert db.execute('SELECT COUNT(*) FROM pairs WHERE kind=?', (kind,)).fetchone()[0] == plan['expected'][kind + '_pairs']
    db.execute('CREATE TABLE triads(id TEXT PRIMARY KEY,payload TEXT,seen INT DEFAULT 0,ready INT DEFAULT 0)')
    previous = ''
    with (out / 'ordered_model_triads.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line); assert previous < row['triad_id']; previous = row['triad_id']
            db.execute('INSERT INTO triads(id,payload) VALUES (?,?)', (row['triad_id'], json.dumps(row)))
    db.execute('CREATE TABLE contexts(guide TEXT,ordinal INT,family TEXT,node TEXT,lo TEXT,hi TEXT,queued TEXT,ab TEXT,PRIMARY KEY(guide,ordinal),UNIQUE(guide,family,node,lo,hi))')
    db.execute('CREATE TABLE context_designs(guide TEXT,ordinal INT,design TEXT,parent INT,empty INT,PRIMARY KEY(guide,ordinal,design))')
    db.execute('CREATE TABLE refs(guide TEXT,ordinal INT,design TEXT,gene TEXT,parent INT,modeled INT,versioned INT,ids INT,pairs INT,ready INT,own INT,both INT,lex INT,lex_both INT,PRIMARY KEY(guide,ordinal,design,gene))')
    contexts = ties = ledger_links = 0
    with source.open() as original, gzip.open(out / 'context_triad_design.jsonl.gz', 'rt') as output:
        for left, right in zip_longest(original, output):
            assert left is not None and right is not None
            work = json.loads(left); actual = json.loads(right)
            assert set(actual) == {'source_design', 'duplicate_pair_design', 'triad_designs'} and actual['source_design'] == work
            context = work['native_context']; guide = context['source_guide']; s = context['source']; parent = context['parent_context_eligible']
            a, b = [(work['duplicate_models'][side]['model_id'], int(work['duplicate_models'][side]['version'])) for side in ['a', 'b']]
            ab = sql_edge(db, a, b, ['a', 'b'], 'primary'); assert actual['duplicate_pair_design'] == ab and ab['pair_key'] == work['duplicate_pair_key']
            queued = work['duplicate_comparison_status']
            if queued == 'queued_distinct_models': assert ab['measurement_disposition'] == 'in_full_primary_measurement_design'
            db.execute('INSERT INTO contexts VALUES (?,?,?,?,?,?,?,?)', (guide, context['source_row_number'], s['family'], s['gene_node'], *sorted([s['gene_a'], s['gene_b']]), queued, ab['measurement_disposition']))
            assert set(actual['triad_designs']) == {'availability', 'sequence_first'}
            for design in ['availability', 'sequence_first']:
                refs, values = work['measurement_designs'][design], actual['triad_designs'][design]; assert len(refs) == len(values)
                db.execute('INSERT INTO context_designs VALUES (?,?,?,?,?)', (guide, context['source_row_number'], design, int(parent), int(not refs)))
                for ref, value in zip(refs, values):
                    reference = ref['reference']; model = reference['reference_model']; triad = None; key = None; reasons = []; versioned = ids = 0; all_pairs = False
                    if model:
                        r = model, int(reference['reference_version']); models = [a, b, r]; key = hash_value(models)
                        edges = dict(ab=ab, ar=sql_edge(db, r, a, ['reference', 'a'], 'reference'), br=sql_edge(db, r, b, ['reference', 'b'], 'reference'))
                        versioned = len(set(models)); ids = len(set(m[0] for m in models))
                        if a == b: reasons.append('duplicate_model_identity')
                        if a == r: reasons.append('a_reference_model_identity')
                        if b == r: reasons.append('b_reference_model_identity')
                        if ids < 3: reasons.append('shared_model_id_across_roles')
                        for label in ['ab', 'ar', 'br']:
                            if edges[label]['measurement_disposition'].startswith('outside_'): reasons.append(label + '_' + edges[label]['measurement_disposition'])
                        all_pairs = edges['ab']['measurement_disposition'] == 'in_full_primary_measurement_design' and edges['ar']['measurement_disposition'] == edges['br']['measurement_disposition'] == 'in_full_reference_measurement_design'
                        triad = dict(triad_id=key, role_order=['a', 'b', 'reference'], models=[list(m) for m in models], distinct_versioned_models=versioned,
                                     distinct_model_ids=ids, edges=edges, all_three_pairs_in_measurement_design=all_pairs, model_work_exclusion_reasons=sorted(reasons))
                        found = db.execute('SELECT payload FROM triads WHERE id=?', (key,)).fetchone(); assert found
                        exported_triad = json.loads(found[0]); assert set(exported_triad) == set(triad) | {'logical_reference_links', 'source_design_ready_links'}
                        assert all(exported_triad[k] == v for k, v in triad.items())
                        for side, label, focal in [('a', 'ar', a), ('b', 'br', b)]:
                            side_work = ref['side_work'][side]; e = edges[label]
                            assert side_work['pair_key'] == e['pair_key'] and side_work['measurement_disposition'] == e['measurement_disposition']
                            if e['current_pair_endpoints'] is not None: assert e['current_pair_endpoints'][side_work['current_pair_focal_endpoint']] == list(focal)
                    else:
                        assert reference['reference_version'] == ''; reasons.append('reference_model_missing')
                    if not parent: reasons.append('parent_context_excluded')
                    if queued != 'queued_distinct_models': reasons.append('duplicate_comparison_' + queued)
                    ready = bool(model and versioned == ids == 3 and all_pairs and parent and queued == 'queued_distinct_models')
                    own = ready and reference['native_coorthology'][guide] == 'both'
                    both = ready and reference['native_coorthology']['profile'] == reference['native_coorthology']['mafft'] == 'both'
                    expected = dict(reference_gene=reference['reference_gene'], triad_id=key, three_distinct_versioned_models=versioned == 3,
                                    three_distinct_model_ids=ids == 3, all_three_pairs_in_measurement_design=bool(all_pairs),
                                    source_design_ready_for_correspondence=ready, source_design_ready_native_own=own,
                                    source_design_ready_native_both_guides=both, work_exclusion_reasons=sorted(set(reasons)))
                    assert value == expected
                    if model: db.execute('UPDATE triads SET seen=seen+1,ready=ready+? WHERE id=?', (int(ready), key))
                    lex = reference['lexical_choice']
                    db.execute('INSERT INTO refs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (guide, context['source_row_number'], design, reference['reference_gene'], int(parent), int(bool(model)), int(versioned == 3), int(ids == 3), int(all_pairs), int(ready), int(own), int(both), int(lex), int(lex and both)))
                    ledger_links += sum(int(x['availability_ledger_link']) for x in ref['side_work'].values()); ties += 1
            contexts += 1
            if contexts % 25000 == 0: print('Independent full SQL three-edge triad design', contexts, '/', plan['expected']['target_contexts'], flush=True)
    guide_counts = dict(db.execute('SELECT guide,COUNT(*) FROM contexts GROUP BY guide'))
    assert contexts == plan['expected']['target_contexts'] and guide_counts == upstream['guide_contexts']
    for guide, n in guide_counts.items(): assert db.execute('SELECT MIN(ordinal),MAX(ordinal),COUNT(*) FROM contexts WHERE guide=?', (guide,)).fetchone() == (1, n, n)
    assert db.execute('SELECT COUNT(*) FROM context_designs').fetchone()[0] == 2 * contexts
    assert ties == plan['expected']['reference_tie_records'] and 2 * ties == plan['expected']['duplicate_reference_links'] and ledger_links == plan['expected']['availability_side_links']
    unique = work_count = 0; model_counts = Counter()
    for payload, seen, ready in db.execute('SELECT payload,seen,ready FROM triads'):
        row = json.loads(payload); assert seen > 0 and row['logical_reference_links'] == seen and row['source_design_ready_links'] == ready
        unique += 1; work_count += int(ready > 0); model_counts[f'versioned={row["distinct_versioned_models"]}|ids={row["distinct_model_ids"]}'] += 1
    duplicates = {f'{g}|{q}|{s}': n for g, q, s, n in db.execute('SELECT guide,queued,ab,COUNT(*) FROM contexts GROUP BY guide,queued,ab')}
    empty = {f'{g}|{d}|parent={p}': n for g, d, p, n in db.execute('SELECT guide,design,parent,COUNT(*) FROM context_designs WHERE empty=1 GROUP BY guide,design,parent')}
    group_fields = ['reference_tie_records', 'modeled_reference_links', 'unmodeled_reference_links', 'three_distinct_versioned_models', 'three_distinct_model_ids', 'all_three_pairs_in_measurement_design', 'source_design_ready_for_correspondence', 'source_design_ready_native_own', 'source_design_ready_native_both_guides', 'lexical_reference_links', 'lexical_source_design_ready_native_both_guides']
    groups = {f'{row[0]}|{row[1]}|parent={row[2]}': dict(zip(group_fields, row[3:])) for row in db.execute('SELECT guide,design,parent,COUNT(*),SUM(modeled),SUM(1-modeled),SUM(versioned),SUM(ids),SUM(pairs),SUM(ready),SUM(own),SUM(both),SUM(lex),SUM(lex_both) FROM refs GROUP BY guide,design,parent')}
    summary = dict(target_contexts=contexts, context_design_records=2 * contexts, reference_tie_records=ties, duplicate_reference_links=2 * ties,
                   availability_side_links_preserved=ledger_links, guide_contexts=guide_counts, unique_ordered_model_triads=unique,
                   correspondence_work_triads=work_count, potential_correspondence_states=16 * work_count, triad_model_identity_counts=dict(model_counts),
                   duplicate_design_counts=duplicates, reference_link_counts=groups, empty_reference_context_counts=empty)
    assert set(summary) == set(SUMMARY_FIELDS) and all(receipt[k] == v for k, v in summary.items()); verify(bindings)
    result = dict(status='passed_full_reference_triad_work_design_sql_readback', plan_sha256=sha(args.plan), producer_receipt_sha256=sha(rp),
                  checker_sha256=sha(__file__), **summary, source_hashes=bindings, scientific_eligibility=False,
                  scope='All full context/native/work fields retained unchanged. Original primary/reference SQL endpoint catalogs independently reconstruct all three desired directions, pair keys, exclusions, model identity and source-only eligibility, including all ties/missing models/excluded parents/empty contexts. Every physical triad exhausted and logical link totals checked independently with SQL aggregates. Shares source I/O and summary field names only; no producer projection functions. No native coverage, residue mapping, biological orthology or duplication effect acceptance.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
