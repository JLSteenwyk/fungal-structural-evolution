#!/usr/bin/env python3
"""Reconstruct all full-reference side states and context screen policies with SQL."""
import argparse
import csv
import gzip
import json
import sqlite3
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from reference_context_coverage_sources import load_sources, POLICIES
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); assert not args.output.exists(); plan = json.loads(args.plan.read_text())
    contexts, table, upstream, _, bindings = load_sources(plan); bind(bindings, args.plan)
    root = Path(plan['output']); receipt_path = root / 'receipt.json'; r = json.loads(receipt_path.read_text())
    assert r['status'] == 'complete_full_reference_context_coverage_pending_independent_readback' and r['plan_sha256'] == sha(args.plan)
    assert r['context_policy_flag_order'] == POLICIES and r['screens'] == plan['screens'] and r['source_hashes'] == bindings
    assert set(r['artifacts']) == {'context_coverage.jsonl.gz', 'context_screen_counts.tsv'}
    bind(bindings, receipt_path)
    for name, h in r['artifacts'].items(): bind(bindings, root / name, h)
    verify(bindings)
    db = sqlite3.connect(':memory:')
    db.execute('CREATE TABLE pairs(key TEXT,mask TEXT,payload TEXT,PRIMARY KEY(key,mask))')
    db.execute('CREATE TABLE contexts(guide TEXT,ordinal INT,design TEXT,parent INT,lex_gene INT,lex_model INT,PRIMARY KEY(guide,ordinal,design))')
    db.execute('CREATE TABLE refs(design TEXT,mask TEXT,screen TEXT,lex INT,a INT,b INT,own INT,both INT)')
    with table.open() as handle:
        db.executemany('INSERT INTO pairs VALUES(?,?,?)', ((x['pair_key'], x['mask'], json.dumps(x)) for x in csv.DictReader(handle, delimiter='\t')))
    assert db.execute('SELECT COUNT(*) FROM pairs').fetchone()[0] == 2 * plan['expected']['model_pairs']
    assert db.execute("SELECT COUNT(*) FROM (SELECT key FROM pairs GROUP BY key HAVING COUNT(*)!=2 OR MIN(mask)!='full' OR MAX(mask)!='plddt70')").fetchone()[0] == 0
    guides, policy_counts, side_counts = Counter(), Counter(), Counter(); ties = logical = n = 0
    with contexts.open() as original, gzip.open(root / 'context_coverage.jsonl.gz', 'rt') as output:
        for left, right in zip_longest(original, output):
            assert left is not None and right is not None
            source, actual = json.loads(left), json.loads(right); context = source['native_context']; guide = context['source_guide']; parent = context['parent_context_eligible']
            assert set(actual) == {'source_design', 'coverage_designs'} and actual['source_design'] == source
            assert set(actual['coverage_designs']) == {'availability', 'sequence_first'}
            pending = []; db.execute('DELETE FROM refs')
            for design in ['availability', 'sequence_first']:
                refs = source['measurement_designs'][design]; exported = actual['coverage_designs'][design]
                assert set(exported) == {'references', 'context_policy_flags'} and len(exported['references']) == len(refs)
                lexical = refs[0]['reference'] if refs else None
                db.execute('INSERT INTO contexts VALUES(?,?,?,?,?,?)', (guide, context['source_row_number'], design, int(parent), int(bool(lexical)), int(bool(lexical and lexical['reference_model']))))
                assert set(exported['context_policy_flags']) == {'full', 'plddt70'}
                assert all(set(s) == {spec['id'] for spec in plan['screens']} for s in exported['context_policy_flags'].values())
                for ix, (ref, coverage) in enumerate(zip(refs, exported['references'])):
                    assert ref['reference']['lexical_choice'] == (ix == 0)
                    assert set(coverage) == {'reference_gene', 'side_masks'} and coverage['reference_gene'] == ref['reference']['reference_gene']
                    assert set(coverage['side_masks']) == {'a', 'b'}
                    passed_sides = {}
                    for side in ['a', 'b']:
                        work = ref['side_work'][side]; disposition = work['measurement_disposition']
                        assert set(coverage['side_masks'][side]) == {'full', 'plddt70'}
                        for mask in ['full', 'plddt70']:
                            row = None
                            if disposition == 'in_full_reference_measurement_design':
                                found = db.execute('SELECT payload FROM pairs WHERE key=? AND mask=?', (work['pair_key'], mask)).fetchone(); assert found; row = json.loads(found[0])
                                ends = [(row['model_a'], int(row['version_a'])), (row['model_b'], int(row['version_b']))]
                                assert ends[work['current_pair_focal_endpoint']] == (work['duplicate_model'], work['duplicate_version'])
                                assert ends[1 - work['current_pair_focal_endpoint']] == (work['reference_model'], work['reference_version'])
                            else: assert disposition in ['reference_model_missing', 'identical_model_not_independent', 'outside_full_reference_measurement_design']
                            expected = dict(source_pair_mask_key=[work['pair_key'], mask] if row else None, order_statuses=[], order_native_statuses=[], order_numerical_exclusions=[], screen_pass={}, screen_exclusions={})
                            for order in [0, 1]:
                                expected['order_statuses'].append(row[f'order{order}_status'] if row else None)
                                expected['order_native_statuses'].append(row[f'order{order}_native_status'] if row else None)
                                text = row[f'order{order}_numerical_exclusion_reasons'] if row else None
                                expected['order_numerical_exclusions'].append(text.split(';') if text else [] if row else None)
                            for spec in plan['screens']:
                                name = spec['id']; reasons = [disposition] if row is None else row[name + '_exclusions'].split(';') if row[name + '_exclusions'] else []
                                pass_value = False if row is None else row[name + '_pass'] == '1'
                                if row: assert row[name + '_pass'] in ['0', '1']
                                assert pass_value == (len(reasons) == 0)
                                if pass_value: assert expected['order_statuses'] == ['aligned', 'aligned']
                                expected['screen_pass'][name] = pass_value; expected['screen_exclusions'][name] = reasons
                                passed_sides[side, mask, name] = int(pass_value)
                            observed = coverage['side_masks'][side][mask]
                            assert observed == expected and all(isinstance(x, bool) for x in observed['screen_pass'].values())
                            side_counts[f'{guide}|{design}|parent={int(parent)}|{mask}|{disposition}'] += 1
                    for mask in ['full', 'plddt70']:
                        for spec in plan['screens']:
                            name = spec['id']; native = ref['reference']['native_coorthology']
                            pending.append((design, mask, name, int(ix == 0), passed_sides['a', mask, name], passed_sides['b', mask, name], int(native[guide] == 'both'), int(native['profile'] == native['mafft'] == 'both')))
                    ties += 1; logical += 2
            db.executemany('INSERT INTO refs VALUES(?,?,?,?,?,?,?,?)', pending)
            aggregate = {}
            for design, mask, screen, count, lexical, own, both, both_ties in db.execute(
                    'SELECT design,mask,screen,COUNT(*),SUM(lex AND a AND b),SUM(lex AND a AND b AND own),SUM(lex AND a AND b AND both),SUM(a AND b AND both) FROM refs GROUP BY design,mask,screen'):
                aggregate[design, mask, screen] = [bool(parent and lexical), bool(parent and own), bool(parent and both), bool(parent and both_ties > 0), bool(parent and count > 0 and both_ties == count)]
            for design in ['availability', 'sequence_first']:
                for mask in ['full', 'plddt70']:
                    for spec in plan['screens']:
                        name = spec['id']; expected = aggregate.get((design, mask, name), [False] * len(POLICIES)); actual_flags = actual['coverage_designs'][design]['context_policy_flags'][mask][name]
                        assert actual_flags == expected and all(isinstance(v, bool) for v in actual_flags)
                        for policy, include in zip(POLICIES, expected): policy_counts[guide, design, mask, name, policy] += include
            guides[guide] += 1; n += 1
            if n % 25000 == 0: print('Independent SQL full context coverage', n, '/', plan['expected']['target_contexts'], flush=True)
    assert n == plan['expected']['target_contexts'] and dict(guides) == upstream['guide_contexts']
    assert db.execute('SELECT COUNT(*) FROM contexts').fetchone()[0] == 2 * n
    for guide, count in guides.items():
        for design in ['availability', 'sequence_first']:
            assert db.execute('SELECT MIN(ordinal),MAX(ordinal),COUNT(*) FROM contexts WHERE guide=? AND design=?', (guide, design)).fetchone() == (1, count, count)
    baseline = {(g, d): (count, eligible, genes, models) for g, d, count, eligible, genes, models in db.execute(
        'SELECT guide,design,COUNT(*),SUM(parent),SUM(parent AND lex_gene),SUM(parent AND lex_model) FROM contexts GROUP BY guide,design')}
    expected_rows = []
    for key, count in sorted(policy_counts.items()):
        base = baseline[key[0], key[1]]
        expected_rows.append(dict(zip(['guide', 'design', 'mask', 'screen', 'policy', 'source_contexts', 'parent_eligible_contexts', 'parent_eligible_lexical_gene_present', 'parent_eligible_lexical_model_present', 'passed_contexts'], map(str, [*key, *base, count]))))
    with (root / 'context_screen_counts.tsv').open() as handle: assert list(csv.DictReader(handle, delimiter='\t')) == expected_rows
    assert ties == plan['expected']['reference_tie_records'] and logical == plan['expected']['duplicate_reference_links'] and sum(side_counts.values()) == 2 * logical
    for cell, count in upstream['measurement_disposition_counts'].items():
        guide, design, parent, status = cell.split('|')
        if status != 'no_reference_gene':
            for mask in ['full', 'plddt70']: assert side_counts[f'{guide}|{design}|{parent}|{mask}|{status}'] == count
    summary = dict(target_contexts=n, context_design_records=2 * n, context_design_mask_records=4 * n, reference_tie_records=ties, duplicate_reference_links=logical,
                   availability_side_links=upstream['availability_side_links_checked'], context_screen_rows=24 * n, context_policy_decisions=24 * n * len(POLICIES),
                   side_screen_decisions=12 * logical, summary_rows=len(expected_rows), guide_contexts=dict(guides), side_mask_work_counts=dict(side_counts))
    assert all(r[k] == v for k, v in summary.items()); verify(bindings)
    result = dict(status='passed_full_reference_context_coverage_sql_readback', plan_sha256=sha(args.plan), producer_receipt_sha256=sha(receipt_path), checker_sha256=sha(__file__),
                  **summary, context_policy_flag_order=POLICIES, source_hashes=bindings, scientific_eligibility=False,
                  scope='Every unchanged full source context/design/tie, model endpoint, measured/unmeasured order status/numerical flag, screen pass and exclusion independently reconstructed through SQL pair lookups. All lexical/any/all-tie context policies independently derived with SQL boolean/count aggregates; original parent eligibility gates even measured overlaps. Empty all-tie sets are excluded. Both guide and complete ordinal/design/mask/screen universes, all source work counts and every summary cell/denominator checked. Shares source/proof I/O only, not producer projection or policy logic. No common-triad/prediction/biological orthology/phylogeny/calibration qualification.')
    with args.output.open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2), flush=True)


if __name__ == '__main__': main()
