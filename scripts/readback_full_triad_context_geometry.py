#!/usr/bin/env python3
"""Independently join full context/model roles and reconstruct geometry policies."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from full_triad_context_geometry_sources import load_sources, DEFINITIONS, MASKS, GATES, POLICIES
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def run(plan_path, output):
    plan_path, output = Path(plan_path), Path(output); assert not output.exists()
    plan = json.loads(plan_path.read_text()); contexts, measurements, upstream, bindings = load_sources(plan, plan_path)
    root = Path(plan['output']); rp = root / 'receipt.json'; r = json.loads(rp.read_text())
    assert r['status'] == 'complete_full_triad_context_geometry_pending_independent_readback' and r['plan_sha256'] == sha(plan_path)
    assert r['masks'] == MASKS and r['definitions'] == DEFINITIONS and r['source_gate_order'] == GATES and r['context_policy_order'] == POLICIES
    assert r['screens'] == plan['screens'] and r['source_hashes'] == bindings and r['scope'] == plan['scope'] and r['scientific_eligibility'] is False
    assert set(r['artifacts']) == {'context_geometry.jsonl.gz', 'context_geometry_counts.tsv'}
    bind(bindings, rp)
    for name, digest in r['artifacts'].items(): bind(bindings, root / name, digest)
    verify(bindings)
    db = sqlite3.connect(':memory:'); screens = [s['id'] for s in plan['screens']]
    db.execute('CREATE TABLE physical(t TEXT,m TEXT,d TEXT,models TEXT,' + ','.join('b' + str(i) + ' INT' for i in range(len(screens))) + ',PRIMARY KEY(t,m,d))')
    db.execute('CREATE TABLE contexts(g TEXT,n INT,d TEXT,parent INT,lex_gene INT,lex_model INT,PRIMARY KEY(g,n,d))')
    db.execute('CREATE TABLE refs(d TEXT,m TEXT,core TEXT,s TEXT,lex INT,ready INT,own INT,both INT,bits INT)')
    groups = 0
    with gzip.open(measurements, 'rt') as handle:
        batch = []
        for line in handle:
            row = json.loads(line); assert row['mask'] in MASKS[:2] and row['mapping_definition'] in DEFINITIONS
            bits = [row['screens'][sid]['order_pass_bits'] for sid in screens]; assert all(type(v) is int and 0 <= v <= 255 for v in bits)
            batch.append([row['triad_id'], row['mask'], row['mapping_definition'], json.dumps(row['models']), *bits]); groups += 1
            if len(batch) == 10000:
                db.executemany('INSERT INTO physical VALUES(' + ','.join(['?'] * (4 + len(screens))) + ')', batch); batch.clear()
        if batch: db.executemany('INSERT INTO physical VALUES(' + ','.join(['?'] * (4 + len(screens))) + ')', batch)
    measured = db.execute('SELECT COUNT(DISTINCT t) FROM physical').fetchone()[0]
    assert measured == plan['expected']['measured_triads'] and groups == plan['expected']['measured_robustness_groups']
    assert db.execute('SELECT COUNT(*) FROM (SELECT t FROM physical GROUP BY t HAVING COUNT(*)!=4 OR COUNT(DISTINCT m)!=2 OR COUNT(DISTINCT d)!=2 OR COUNT(DISTINCT models)!=1)').fetchone()[0] == 0
    counts, guides, presence = Counter(), Counter(), Counter(); n = ties = 0
    with gzip.open(contexts, 'rt') as original, gzip.open(root / 'context_geometry.jsonl.gz', 'rt') as exported:
        for left, right in zip_longest(original, exported):
            assert left is not None and right is not None
            source, actual = json.loads(left), json.loads(right)
            assert set(actual) == {'source_work_design', 'structural_designs'} and actual['source_work_design'] == source
            assert set(actual['structural_designs']) == {'availability', 'sequence_first'}
            work = source['source_design']; native = work['native_context']; guide, parent = native['source_guide'], native['parent_context_eligible']
            db.execute('DELETE FROM refs'); pending = []
            for design in ['availability', 'sequence_first']:
                output_design = actual['structural_designs'][design]
                assert set(output_design) == {'references', 'context_policy_flags'}
                refs = work['measurement_designs'][design]; designs = source['triad_designs'][design]
                assert len(refs) == len(designs) == len(output_design['references'])
                lexical = refs[0]['reference'] if refs else None
                db.execute('INSERT INTO contexts VALUES(?,?,?,?,?,?)', (guide, native['source_row_number'], design, int(parent), int(bool(lexical)), int(bool(lexical and lexical['reference_model']))))
                for ix, (ref, disposition, observed) in enumerate(zip(refs, designs, output_design['references'])):
                    chosen = ref['reference']; assert chosen['lexical_choice'] == (ix == 0) and disposition['reference_gene'] == chosen['reference_gene']
                    models = [[work['duplicate_models'][role]['model_id'], work['duplicate_models'][role]['version']] for role in ['a', 'b']]
                    if chosen['reference_model']:
                        models.append([chosen['reference_model'], int(chosen['reference_version'])])
                        key = hashlib.sha256(json.dumps(models, separators=(',', ':')).encode()).hexdigest()
                        distinct_versioned = len({tuple(m) for m in models}) == 3; distinct_ids = len({m[0] for m in models}) == 3
                        all_pairs = source['duplicate_pair_design']['measurement_disposition'] == 'in_full_primary_measurement_design' and all(ref['side_work'][side]['measurement_disposition'] == 'in_full_reference_measurement_design' for side in ['a', 'b'])
                    else: key = None; distinct_versioned = distinct_ids = all_pairs = False
                    assert disposition['triad_id'] == key and disposition['three_distinct_versioned_models'] == distinct_versioned and disposition['three_distinct_model_ids'] == distinct_ids and disposition['all_three_pairs_in_measurement_design'] == all_pairs
                    ready = bool(key and distinct_versioned and distinct_ids and all_pairs and parent and work['duplicate_comparison_status'] == 'queued_distinct_models')
                    own = bool(ready and chosen['native_coorthology'][guide] == 'both'); both = bool(ready and chosen['native_coorthology']['profile'] == chosen['native_coorthology']['mafft'] == 'both')
                    gates = [ready, own, both]; assert [disposition[gate] for gate in GATES] == gates
                    physical = {(m, d): (json.loads(payload), bits) for m, d, payload, *bits in db.execute('SELECT m,d,models,' + ','.join('b' + str(i) for i in range(len(screens))) + ' FROM physical WHERE t=?', (key,))}
                    present = bool(physical); assert not physical or len(physical) == 4
                    if present: assert all(payload == models for payload, _ in physical.values())
                    if ready: assert present
                    order_bits = {}; qualified = {}
                    for mask in MASKS:
                        order_bits[mask] = {}; qualified[mask] = {}
                        for definition in DEFINITIONS:
                            cells = {}; q = {}
                            for i, sid in enumerate(screens):
                                bits = None
                                if present:
                                    bits = physical[mask, definition][1][i] if mask != 'both_masks' else physical['full', definition][1][i] & physical['plddt70', definition][1][i]
                                cells[sid] = bits; q[sid] = [bool(bits == 255 and gate) for gate in gates]
                                pending.append((design, mask, definition, sid, int(ix == 0), *map(int, gates), bits))
                            order_bits[mask][definition] = cells; qualified[mask][definition] = q
                    expected = dict(reference_gene=chosen['reference_gene'], triad_id=key, physical_measurement_present=present,
                                    physical_robustness_keys=[[key, mask, definition] for mask in MASKS[:2] for definition in DEFINITIONS] if present else None,
                                    source_eligibility_flags=gates, order_pass_bits=order_bits, qualified_all_order_flags=qualified)
                    assert observed == expected and type(observed['physical_measurement_present']) is bool
                    assert all(type(v) is bool for v in observed['source_eligibility_flags'])
                    assert all(type(v) is bool for m in MASKS for d in DEFINITIONS for sid in screens for v in observed['qualified_all_order_flags'][m][d][sid])
                    presence[f'{guide}|{design}|parent={int(parent)}|measured={int(present)}|source_ready={int(ready)}'] += 1; ties += 1
            db.executemany('INSERT INTO refs VALUES(?,?,?,?,?,?,?,?,?)', pending)
            aggregates = {}
            for design, mask, definition, sid, count, lexical, own, both, passed in db.execute(
                    'SELECT d,m,core,s,COUNT(*),SUM(lex AND ready AND COALESCE(bits,0)=255),SUM(lex AND own AND COALESCE(bits,0)=255),SUM(lex AND both AND COALESCE(bits,0)=255),SUM(both AND COALESCE(bits,0)=255) FROM refs GROUP BY d,m,core,s'):
                aggregates[design, mask, definition, sid] = [bool(lexical), bool(own), bool(both), passed > 0, count > 0 and passed == count]
            for design in ['availability', 'sequence_first']:
                expected_policies = {}
                for mask in MASKS:
                    expected_policies[mask] = {}
                    for definition in DEFINITIONS:
                        expected_policies[mask][definition] = {}
                        for sid in screens:
                            flags = aggregates.get((design, mask, definition, sid), [False] * 5)
                            expected_policies[mask][definition][sid] = flags
                            for policy, flag in zip(POLICIES, flags): counts[guide, design, mask, definition, sid, policy] += flag
                observed = actual['structural_designs'][design]['context_policy_flags']; assert observed == expected_policies
                assert all(type(v) is bool for m in MASKS for d in DEFINITIONS for sid in screens for v in observed[m][d][sid])
            n += 1; guides[guide] += 1
            if n % 25000 == 0: print('Independent SQL full context geometry', n, '/', plan['expected']['target_contexts'], flush=True)
    assert n == plan['expected']['target_contexts'] and ties == plan['expected']['reference_tie_records'] and dict(guides) == upstream['guide_contexts']
    assert db.execute('SELECT COUNT(*) FROM contexts').fetchone()[0] == 2 * n
    for guide, count in guides.items():
        for design in ['availability', 'sequence_first']:
            assert db.execute('SELECT MIN(n),MAX(n),COUNT(*) FROM contexts WHERE g=? AND d=?', (guide, design)).fetchone() == (1, count, count)
    baseline = {(g, d): (count, parent, genes, models) for g, d, count, parent, genes, models in db.execute('SELECT g,d,COUNT(*),SUM(parent),SUM(parent AND lex_gene),SUM(parent AND lex_model) FROM contexts GROUP BY g,d')}
    fields = ['guide', 'design', 'mask', 'mapping_definition', 'screen', 'policy', 'source_contexts', 'parent_eligible_contexts', 'parent_eligible_lexical_gene_present', 'parent_eligible_lexical_model_present', 'passed_contexts']
    expected_rows = [dict(zip(fields, map(str, [*key, *baseline[key[:2]], count]))) for key, count in sorted(counts.items())]
    with (root / 'context_geometry_counts.tsv').open() as handle: assert list(csv.DictReader(handle, delimiter='\t')) == expected_rows
    summary = dict(target_contexts=n, context_design_records=2 * n, reference_tie_records=ties, duplicate_reference_links=2 * ties,
                   logical_reference_screen_decisions=ties * len(MASKS) * len(DEFINITIONS) * len(screens), context_screen_states=n * 2 * len(MASKS) * len(DEFINITIONS) * len(screens),
                   context_policy_decisions=n * 2 * len(MASKS) * len(DEFINITIONS) * len(screens) * len(POLICIES), summary_rows=len(counts), guide_contexts=dict(guides),
                   reference_measurement_presence_counts=dict(presence), measured_triads=measured, measured_robustness_groups=groups)
    assert summary['duplicate_reference_links'] == plan['expected']['duplicate_reference_links'] and all(r[k] == value for k, value in summary.items()); verify(bindings)
    result = dict(status='passed_full_triad_context_geometry_sql_readback', plan_sha256=sha(plan_path), producer_receipt_sha256=sha(rp), checker_sha256=sha(__file__),
                  **summary, masks=MASKS, definitions=DEFINITIONS, source_gate_order=GATES, context_policy_order=POLICIES, source_hashes=bindings, scientific_eligibility=False,
                  scope='All unchanged source contexts/designs/ties joined independently by original ordered model roles and SHA keys. Recomputes source readiness from original parent, model identities, all3edge design and native guides. All nullable physical8order bitmaps, both-mask intersections and gated all-order reference flags checked. Independently reconstructs lexical/any/all-tie policy vectors and every summary/denominator with SQL; empty ties excluded. Physical measured overlaps cannot override source parent/model/orthology gates. Shares source/proof I/O only, no producer linkage/policy code; not calibrated evolutionary inference.')
    with output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', required=True, type=Path); parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(); run(args.plan, args.output)
