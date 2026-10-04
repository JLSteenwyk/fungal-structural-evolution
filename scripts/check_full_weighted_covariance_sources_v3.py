#!/usr/bin/env python3
"""Complete synthetic source/grid checks with semantically rehashed corruptions."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import sqlite3

import numpy as np
import pyarrow.parquet as pq
from scipy import sparse

from ancestral_chain_attempt import sha
from check_full_covariance_qualification import setup as original_setup
from check_full_inverse_reuse_weights import closure, write
from covariance_exact_folds_v2 import certificate, variance_map
from full_expanded_model_design_sources import OUTCOMES, array_digest, digest
from full_weighted_covariance_sources_v2 import MODES, STATUSES, basis, cohort_controls, load
from inverse_reuse_weight_controls import POLICIES, calculate
from nonuniform_covariance_cone import record as cone_record
from reference_measurement_union_sources import verify
from run_full_inverse_reuse_weights import record as weight_record, RECIPE
from run_full_weighted_covariance_source_census_v2 import run


def rejected(action):
    try: action()
    except (AssertionError, ValueError, KeyError, FileExistsError, StopIteration, sqlite3.IntegrityError): return
    raise AssertionError('Malformed complete weighted source accepted')


def setup(root):
    root.mkdir(); oldp = original_setup(root); old = json.loads(oldp.read_text())
    configs = {k: json.loads(Path(old[k + '_plan']).read_text()) for k in ['operator', 'design', 'inputs']}
    roots = {k: Path(v['output']) for k, v in configs.items()}
    bank = configs['operator']; covp = Path(bank['covariance_plan']); covconfig = json.loads(covp.read_text())
    roots['covariance'] = Path(covconfig['output'])
    cohorts = [c for c in json.loads((roots['design'] / 'cohort_manifest.json').read_text()) if c['records'] > 0]
    write(roots['design'] / 'cohort_manifest.json', cohorts); cids = {c['cohort_id'] for c in cohorts}
    ds = roots['design'] / 'unique_designs.jsonl'
    designs = [json.loads(l) for l in ds.read_text().splitlines() if json.loads(l)['cohort_id'] in cids]
    ds.write_text(''.join(json.dumps(d) + '\n' for d in designs)); dids = {d['design_id'] for d in designs}
    fp = roots['design'] / 'unique_fit_inputs.jsonl'
    oldfits = [json.loads(l) for l in fp.read_text().splitlines() if json.loads(l)['design_id'] in dids]
    parts = json.loads((roots['inputs'] / 'partition_manifest.json').read_text())
    arrays = {(p['mask'], p['order_contrast']): pq.read_table(roots['inputs'] / p['path']).to_pydict() for p in parts}
    byc = {c['cohort_id']: c for c in cohorts}; byd = {d['design_id']: d for d in designs}
    for f in oldfits:
        d = byd[f['design_id']]; c = byc[d['cohort_id']]
        with np.load(roots['design'] / c['path']) as a: selected = a['case_rows']
        y = np.asarray(arrays[c['mask'], d['order_contrast']][f['outcome']])[selected]
        f.update(cohort_id=c['cohort_id'], records=len(y), trees=old['trees'],
            response_sha256=array_digest(y, '<f8'), response_min=float(y.min()), response_max=float(y.max()))
    fp.write_text(''.join(json.dumps(f) + '\n' for f in oldfits))
    sp = roots['design'] / 'model_settings.tsv.gz'
    lines = gzip.decompress(sp.read_bytes()).decode().splitlines(); fields = lines[0].split('\t'); index = fields.index('cohort_id')
    lines = [lines[0]] + [l for l in lines[1:] if l.split('\t')[index] in cids]
    sp.write_bytes(gzip.compress(('\n'.join(lines) + '\n').encode(), mtime=0))
    n = 24; count = len(cohorts); occurrences = sum(c['records'] for c in cohorts)
    pp = {k: Path(old[k + '_plan']) for k in ['operator', 'design', 'inputs']}
    cp = {k: Path(old[k + '_completion']) for k in ['operator', 'design', 'inputs']}
    pp['covariance'] = covp; cp['covariance'] = Path(bank['covariance_completion'])
    configs['design'].update(inputs_plan=str(pp['inputs']), inputs_completion=str(cp['inputs']))
    write(pp['design'], configs['design'])
    for k in ['exact', 'cone', 'weights']:
        roots[k] = root / k; roots[k].mkdir(); pp[k] = root / (k + '.plan.json'); cp[k] = root / (k + '.completed.json')
    write(pp['exact'], dict(output=str(roots['exact'])))
    write(pp['cone'], dict(output=str(roots['cone']), parent_plan=str(pp['exact']), parent_completion=str(cp['exact'])))
    write(pp['weights'], dict(output=str(roots['weights']),
        **{k + '_plan': str(pp[k]) for k in ['operator', 'design', 'covariance', 'cone']},
        **{k + '_completion': str(cp[k]) for k in ['operator', 'design', 'covariance', 'cone']}))
    operators = {(r['mode'], r['kind']): sparse.load_npz(roots['operator'] / r['path']).tocsr()
        for r in json.loads((roots['operator'] / 'operator_manifest.json').read_text())}
    labels = np.load(roots['operator'] / 'block_labels.npy'); manifest = []; parents = []; cones = []
    wc = digest(['synthetic-weight-contract']); cc = digest(['synthetic-cone-contract'])
    (roots['weights'] / 'arrays').mkdir(); (roots['weights'] / 'cohorts').mkdir()
    for c in cohorts:
        with np.load(roots['design'] / c['path']) as a: selected = a['case_rows']
        for mode in MODES:
            ops = {k: operators[mode, k][selected].tocsr() for k in ['target_node', 'background_node', 'model_pair', 'gene', 'model', 'family']}
            family = operators['contrast', 'family_intercept'][selected].tocsr()
            proof = certificate(ops, family, mode, selected)
            parent = dict(cohort_id=c['cohort_id'], cohort_rows_sha256=c['case_rows_sha256'],
                ordered_case_ids_sha256=c['ordered_case_ids_sha256'], certificate=proof,
                variance_map=variance_map(proof['relations'], mode), raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
            parents.append(parent); cones.append(cone_record(parent, cc))
        keys = [np.asarray([str(i // 3) for i in selected]), np.asarray([str(i // 4) for i in selected]), labels[selected]]
        counts, weights, diagonal, groups = calculate(keys)
        ap = roots['weights'] / 'arrays' / (c['cohort_id'] + '.npz')
        np.savez_compressed(ap, case_rows=selected, reuse_counts=counts, weights=weights, reciprocal_diagonals=diagonal)
        jp = roots['weights'] / 'cohorts' / (c['cohort_id'] + '.json')
        write(jp, weight_record(c, wc, counts, weights, diagonal, groups, selected, ap))
        manifest.append(dict(cohort_id=c['cohort_id'], records=len(selected), record_path=str(jp.relative_to(roots['weights'])),
            record_sha256=sha(jp), array_path=str(ap.relative_to(roots['weights'])), array_sha256=sha(ap)))
    for k, name, records in [('exact', 'cohort_certificates.jsonl', parents), ('cone', 'cohort_cones.jsonl', cones)]:
        (roots[k] / name).write_text(''.join(json.dumps(r) + '\n' for r in records))
    write(roots['cone'] / 'receipt.json', dict(source_contract=cc))
    write(roots['weights'] / 'control_manifest.json', manifest)
    write(roots['weights'] / 'receipt.json', dict(source_contract=wc, recipe=RECIPE, policies=POLICIES))
    summaries = {k: dict(logical_cases=n) for k in STATUSES}
    summaries['design'].update(unique_cohorts=count, unique_designs=len(designs), model_setting_rows=len(lines) - 1, trees=old['trees'])
    for k in ['exact', 'cone', 'weights']: summaries[k].update(cohorts=count)
    for k in ['exact', 'cone']: summaries[k].update(certificates=count * 2)
    summaries['weights'].update(policies=POLICIES, residual_diagonal_prepared=True,
        raw_reml_basis_qualification_complete=False, nonuniform_weighting_accepted=False, case_row_occurrences=occurrences)
    for k in STATUSES:
        paths = [pp[k], *[p for p in roots[k].rglob('*') if p.is_file() and p not in [cp[k]] and not p.name.endswith('archive.json') and p.name != 'completion.json']]
        closure(cp[k], STATUSES[k], paths, summaries[k])
    plan = {k + '_plan': str(pp[k]) for k in STATUSES}; plan.update({k + '_completion': str(cp[k]) for k in STATUSES})
    modules = ['full_weighted_covariance_sources_v2', 'run_full_weighted_covariance_source_census_v2', 'check_full_weighted_covariance_sources_v3',
        'full_covariance_qualification_sources', 'full_exact_covariance_sources', 'full_expanded_model_design_sources',
        'full_expanded_model_input_sources', 'inverse_reuse_weight_controls', 'nonuniform_covariance_cone',
        'reduced_covariance_basis', 'covariance_exact_folds_v2', 'reference_measurement_union_sources', 'run_full_inverse_reuse_weights']
    plan.update(pins={'scripts/' + m + '.py': sha('scripts/' + m + '.py') for m in modules}, trees=old['trees'],
        output=str(root / 'census'), expected=dict(logical_cases=n, cohorts=count, designs=len(designs), fit_inputs=len(oldfits),
        settings=len(lines) - 1, case_row_occurrences=occurrences,
        numerical_audit_rows=len(designs) * 40, setting_audit_links=(len(lines) - 1) * 40),
        scope='Complete declared synthetic five nonempty original cohorts, all thirty designs, two outcomes, two modes, five trees and four controls. Source closures/journals are synthetic. No real numerical qualification, fit or biological pilot.')
    path = root / 'weighted.plan.json'; write(path, plan); return path


def rehash(plan, kind, paths):
    cp = Path(plan[kind + '_completion']); c = json.loads(cp.read_text()); ap = Path(c['full_hash_archive']); a = json.loads(ap.read_text())
    for p in paths: a['source_hashes'][str(p)] = sha(p)
    write(ap, a); c['full_hash_archive_sha256'] = sha(ap); write(cp, c)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True); a = p.parse_args(); assert not a.receipt.exists()
    a.output.mkdir(exist_ok=False); pp = setup(a.output / 'fixture'); plan = json.loads(pp.read_text())
    produced = run(pp); checked = run(pp, True)
    assert (produced['cohorts'], produced['designs'], produced['fit_inputs'], produced['settings'],
        produced['expected_numerical_audit_rows'], produced['expected_setting_audit_links']) == (5, 150, 300, 600, 6000, 24000)
    rejected(lambda: run(pp)); rejected(lambda: run(pp, True))
    source, bindings = load(plan, pp)
    routechecks = 0
    for c, entry in zip(source['cohorts'], source['controls']):
        selected, diagonals, control = cohort_controls(source, c, entry)
        for d in [np.ones(len(selected)), np.full(len(selected), 2.), 1. + np.arange(len(selected)) * 2.**-40]:
            for mode in MODES:
                route = basis(source, c, mode, d)
                assert ('target_node' not in route['names']) == bool(np.array_equal(d, np.ones(len(d))))
                routechecks += 1
    root = Path(plan['output']); rp = root / 'receipt.json'; rb = root / 'readback.json'; data = root / 'cohort_source_census.jsonl.gz'
    saved = {p: p.read_bytes() for p in [rp, rb, data]}; rb.unlink()
    outputcases = ['missing_cohort', 'duplicate_cohort', 'foreign_design', 'missing_fit', 'drop_control', 'drop_target_identity', 'alter_diagonal_hash', 'promote_science']
    for name in outputcases:
        records = [json.loads(l) for l in gzip.decompress(saved[data]).decode().splitlines()]
        if name == 'missing_cohort': records.pop()
        elif name == 'duplicate_cohort': records[-1] = records[0]
        elif name == 'foreign_design': records[0]['design_ids'][0] = 'foreign'
        elif name == 'missing_fit': records[0]['fit_input_ids'].pop()
        elif name == 'drop_control': records[0]['policy_basis_choices'].pop()
        elif name == 'drop_target_identity':
            r = next(r for rec in records for r in rec['policy_basis_choices'] if 'target_node' in r['names']); r['names'].remove('target_node')
        elif name == 'alter_diagonal_hash': records[0]['policy_basis_choices'][0]['diagonal_sha256'] = 'foreign'
        else: records[0]['scientific_eligibility'] = True
        data.write_bytes(gzip.compress((''.join(json.dumps(r) + '\n' for r in records)).encode(), mtime=0))
        receipt = json.loads(saved[rp]); receipt['artifacts'][data.name] = sha(data); write(rp, receipt)
        rejected(lambda: run(pp, True)); rp.write_bytes(saved[rp]); data.write_bytes(saved[data])
    # Rehash source closures deliberately; identity/schema/membership guards
    # must reject semantic changes even when all affected digests agree.
    sourcecases = ['missing_tree', 'foreign_case_order', 'wrong_pattern_row', 'foreign_raw_design', 'wrong_active_columns',
        'wrong_response_hash', 'foreign_setting_link', 'foreign_cone', 'foreign_weight_recipe', 'wrong_control_membership']
    originals = {p: p.read_bytes() for p in a.output.rglob('*') if p.is_file()}
    for name in sourcecases:
        if name == 'missing_tree':
            mutated = dict(plan); mutated['trees'] = plan['trees'][:-1]; write(pp, mutated)
        else:
            kind = 'operator' if name == 'foreign_case_order' else 'covariance' if name == 'wrong_pattern_row' else 'cone' if name == 'foreign_cone' else 'weights' if name.startswith('foreign_weight') or name == 'wrong_control_membership' else 'design'
            sr = source['roots'][kind]
            if name == 'foreign_case_order':
                q = sr / 'case_ids.json'; write(q, list(reversed(json.loads(q.read_text()))))
            elif name == 'wrong_pattern_row':
                q = sr / 'case_covariance_index.tsv.gz'; text = gzip.decompress(q.read_bytes()).decode().splitlines(); fields = text[0].split('\t')
                cells = text[1].split('\t'); cells[fields.index('species_pattern_row')] = '-1'; text[1] = '\t'.join(cells)
                q.write_bytes(gzip.compress(('\n'.join(text) + '\n').encode(), mtime=0))
            elif name in ['foreign_raw_design', 'wrong_active_columns', 'wrong_response_hash', 'foreign_cone']:
                q = sr / ('cohort_cones.jsonl' if name == 'foreign_cone' else 'unique_fit_inputs.jsonl' if name == 'wrong_response_hash' else 'unique_designs.jsonl')
                records = [json.loads(l) for l in q.read_text().splitlines()]
                if name == 'foreign_raw_design': records[0]['raw_design_sha256'] = 'foreign'
                elif name == 'wrong_active_columns': records[0]['active_column_indices'] = list(reversed(records[0]['active_column_indices']))
                elif name == 'wrong_response_hash': records[0]['response_sha256'] = 'foreign'
                else: records[0]['parent_exact_certificate_sha256'] = 'foreign'
                q.write_text(''.join(json.dumps(r) + '\n' for r in records))
            elif name == 'foreign_setting_link':
                q = sr / 'model_settings.tsv.gz'; text = gzip.decompress(q.read_bytes()).decode().splitlines(); fields = text[0].split('\t')
                cells = text[1].split('\t'); cells[fields.index('fit_input_id')] = 'foreign'; text[1] = '\t'.join(cells)
                q.write_bytes(gzip.compress(('\n'.join(text) + '\n').encode(), mtime=0))
            else:
                q = sr / source['controls'][0]['record_path']; r = json.loads(q.read_text())
                if name == 'foreign_weight_recipe': r['recipe']['selection_record_multiplicity_used'] = True
                else: r['cohort_rows_sha256'] = 'foreign'
                write(q, r); mp = sr / 'control_manifest.json'; m = json.loads(mp.read_text()); m[0]['record_sha256'] = sha(q); write(mp, m)
                rehash(plan, kind, [mp])
            rehash(plan, kind, [q])
        mutated = json.loads(pp.read_text()); mutated['output'] = str(a.output / ('rejected-source-' + name))
        mutation_plan = a.output / ('rejected-source-' + name + '.json'); write(mutation_plan, mutated)
        rejected(lambda: run(mutation_plan))
        for q, b in originals.items(): q.write_bytes(b)
    replay = run(pp, True)
    # The timestamp differs; compare all substantive reader evidence.
    previous = json.loads(saved[rb]); current = dict(replay)
    previous.pop('checked_utc'); current.pop('checked_utc'); assert previous == current
    verify(replay['source_hashes'])
    result = dict(status='passed_complete_declared_four_control_source_census_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), synthetic_logical_cases=24, synthetic_cohorts=5,
        synthetic_designs=150, synthetic_fit_inputs=300, synthetic_settings=600,
        prospective_audits=6000, prospective_setting_links=24000, exact_uniform_constant_near_uniform_route_checks=routechecks,
        rehashed_output_cases_rejected=outputcases, rehashed_source_cases_rejected=sourcecases,
        completed_stage_restart_refusals=True, substantive_positive_restoration=True,
        source_and_journal_fixtures_synthetic=True, numerical_audits_computed=0,
        scientific_eligibility=False, source_hashes={**replay['source_hashes'], str(root / 'readback.json'): sha(root / 'readback.json')})
    with a.receipt.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__': main()
