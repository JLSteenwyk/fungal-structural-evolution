"""Fresh, scoped original inputs for all four covariance sensitivity controls.

Closed operator identities authorize the basis. No saved numerical audit or
uniform error envelope is consumed. Group membership itself was independently
verified by the closed SQL/Fraction control stage; this loader checks its actual
saved arrays, recipe, cohort identity and reciprocal arithmetic again.
"""
import csv
import gzip
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy import sparse

from ancestral_chain_attempt import sha
from full_covariance_qualification_sources import design_matrix, jsonl
from full_exact_covariance_sources import closed_subset, rows
from full_expanded_model_design_sources import AXES, OUTCOMES, array_digest, digest
from full_expanded_model_input_sources import MASKS, ORDERS, NUISANCE, schema
from inverse_reuse_weight_controls import ARRAYS, CLAIMS, POLICIES, SCHEMA as WEIGHT_SCHEMA, policy_records
from nonuniform_covariance_cone import record as cone_record
from reduced_covariance_basis import validate_certificate
from reference_measurement_union_sources import bind, verify
from run_full_inverse_reuse_weights import RECIPE

SCHEMA = 'full-four-control-covariance-inputs-v1'
MODES = ['signed', 'unsigned']
STATUSES = dict(operator='complete_verified_full_entity_operator_bank',
    design='complete_verified_full_expanded_model_designs',
    inputs='complete_verified_full_expanded_model_inputs',
    covariance='complete_verified_full_expanded_covariance',
    exact='complete_verified_full_exact_uniform_covariance_folds_v2',
    cone='complete_verified_full_positive_diagonal_covariance_cones_v1',
    weights='complete_verified_full_original_cohort_inverse_reuse_controls_v1')


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path); verify(bindings)
    configs = {k: json.loads(Path(plan[k + '_plan']).read_text()) for k in STATUSES}
    roots = {k: Path(v['output']) for k, v in configs.items()}
    assert configs['operator']['covariance_plan'] == plan['covariance_plan']
    assert configs['operator']['covariance_completion'] == plan['covariance_completion']
    assert configs['design']['inputs_plan'] == plan['inputs_plan']
    assert configs['design']['inputs_completion'] == plan['inputs_completion']
    assert configs['cone']['parent_plan'] == plan['exact_plan']
    assert configs['cone']['parent_completion'] == plan['exact_completion']
    for k in ['operator', 'design', 'covariance', 'cone']:
        assert configs['weights'][k + '_plan'] == plan[k + '_plan']
        assert configs['weights'][k + '_completion'] == plan[k + '_completion']
    cohorts = json.loads((roots['design'] / 'cohort_manifest.json').read_text())
    operators = json.loads((roots['operator'] / 'operator_manifest.json').read_text())
    partitions = json.loads((roots['inputs'] / 'partition_manifest.json').read_text())
    controls = json.loads((roots['weights'] / 'control_manifest.json').read_text())
    paths = dict(
        design=[roots['design'] / n for n in ['receipt.json', 'cohort_manifest.json', 'unique_designs.jsonl',
            'unique_fit_inputs.jsonl', 'model_settings.tsv.gz']] + [roots['design'] / c['path'] for c in cohorts],
        operator=[roots['operator'] / n for n in ['case_ids.json', 'block_labels.npy', 'operator_manifest.json']]
            + [roots['operator'] / r['path'] for r in operators],
        inputs=[roots['inputs'] / 'partition_manifest.json'] + [roots['inputs'] / r['path'] for r in partitions],
        covariance=[roots['covariance'] / 'case_covariance_index.tsv.gz']
            + [roots['covariance'] / (t + '.npz') for t in plan['trees']],
        exact=[roots['exact'] / 'cohort_certificates.jsonl'],
        cone=[roots['cone'] / 'receipt.json', roots['cone'] / 'cohort_cones.jsonl'],
        weights=[roots['weights'] / 'receipt.json', roots['weights'] / 'control_manifest.json']
            + [roots['weights'] / r[k] for r in controls for k in ['record_path', 'array_path']])
    closed = {}
    for k in STATUSES:
        closed[k] = closed_subset(plan[k + '_completion'], STATUSES[k],
            [Path(plan[k + '_plan']), *paths[k]], bindings)
    n = plan['expected']['logical_cases']; count = plan['expected']['cohorts']
    assert all(c['logical_cases'] == n for c in closed.values())
    assert closed['design']['unique_cohorts'] == closed['exact']['cohorts'] == closed['cone']['cohorts'] == closed['weights']['cohorts'] == count
    assert len(cohorts) == len(controls) == count
    assert [c['cohort_id'] for c in cohorts] == sorted({c['cohort_id'] for c in cohorts})
    assert [c['cohort_id'] for c in controls] == [c['cohort_id'] for c in cohorts]
    assert closed['design']['unique_designs'] == plan['expected']['designs'] == count * 30
    assert closed['design']['model_setting_rows'] == plan['expected']['settings']
    assert closed['weights']['policies'] == POLICIES and closed['weights']['residual_diagonal_prepared'] is True
    assert closed['weights']['raw_reml_basis_qualification_complete'] is closed['weights']['nonuniform_weighting_accepted'] is False
    assert closed['weights']['case_row_occurrences'] == plan['expected']['case_row_occurrences']
    assert closed['design']['trees'] == configs['design']['trees'] == configs['operator']['trees'] == plan['trees']
    assert len(plan['trees']) == 5 and len(set(plan['trees'])) == 5
    ids = json.loads((roots['operator'] / 'case_ids.json').read_text())
    labels = np.load(roots['operator'] / 'block_labels.npy', allow_pickle=False)
    assert len(ids) == len(set(ids)) == n and labels.shape == (n,) and labels.dtype.kind in 'US'
    expected = {(m, k) for m in MODES for k in ['target_node', 'background_node', 'model_pair', 'gene', 'model', 'family']} | {('contrast', 'family_intercept')}
    assert len(operators) == 13 and {(r['mode'], r['kind']) for r in operators} == expected
    bank = {}
    for r in operators:
        p = roots['operator'] / r['path']; assert bindings[str(p)] == r['sha256']
        z = sparse.load_npz(p).tocsr()
        assert list(z.shape) == r['shape'] and z.nnz == r['nnz'] and z.shape[0] == n
        assert z.has_canonical_format and np.isfinite(z.data).all()
        bank[r['mode'], r['kind']] = z
    assert [(r['mask'], r['order_contrast']) for r in partitions] == [(m, o) for m in MASKS for o in ORDERS]
    columns = ['case_id', 'input_id', *dict.fromkeys(sum(AXES.values(), []) + NUISANCE + OUTCOMES)]
    arrays = {}
    for r in partitions:
        p = roots['inputs'] / r['path']; assert bindings[str(p)] == r['sha256']
        pf = pq.ParquetFile(p); assert pf.schema_arrow == schema() and pf.metadata.num_rows == n
        t = pf.read(columns=columns); assert t['case_id'].to_pylist() == ids
        arrays[r['mask'], r['order_contrast']] = {k: np.asarray(t[k].to_pylist(), dtype='S64' if k in columns[:2] else float) for k in columns}
    with gzip.open(roots['covariance'] / 'case_covariance_index.tsv.gz', 'rt') as f:
        cases = list(csv.DictReader(f, delimiter='\t'))
    assert [c['case_id'] for c in cases] == ids
    pattern_rows = np.asarray([int(c['species_pattern_row']) for c in cases], dtype=np.int64)
    components = [c['family_component'] for c in cases]
    assert labels.astype(str).tolist() == components
    factors = {}
    for tree in plan['trees']:
        with np.load(roots['covariance'] / (tree + '.npz'), allow_pickle=False) as a:
            f = a['factor']; assert f.ndim == 2 and np.isfinite(f).all()
            assert np.all((pattern_rows >= 0) & (pattern_rows < len(f)))
            factors[tree] = f[pattern_rows]
    parents = {}; cones = {}; by_id = {c['cohort_id']: c for c in cohorts}
    for r in jsonl(roots['exact'] / 'cohort_certificates.jsonl'):
        mode = r['certificate']['loading_mode']; key = r['cohort_id'], mode
        assert key not in parents and mode in MODES
        c = by_id[r['cohort_id']]; validate_certificate(r)
        assert r['certificate']['records'] == c['records']
        assert r['cohort_rows_sha256'] == c['case_rows_sha256']
        assert r['ordered_case_ids_sha256'] == c['ordered_case_ids_sha256']
        parents[key] = r
    cone_receipt = json.loads((roots['cone'] / 'receipt.json').read_text())
    for r in jsonl(roots['cone'] / 'cohort_cones.jsonl'):
        key = r['cohort_id'], r['loading_mode']; assert key not in cones
        assert r == cone_record(parents[key], cone_receipt['source_contract'])
        cones[key] = r
    assert set(parents) == set(cones) == {(c['cohort_id'], m) for c in cohorts for m in MODES}
    assert closed['exact']['certificates'] == closed['cone']['certificates'] == 2 * count
    receipt = json.loads((roots['weights'] / 'receipt.json').read_text())
    assert receipt['recipe'] == RECIPE and receipt['policies'] == POLICIES
    contract = digest(dict(schema=SCHEMA, source_hashes=bindings, policies=POLICIES,
        basis='named uniform certificate only for actual all-one D; positive-diagonal cone otherwise',
        numerical_audits_inherited=False, scope='fresh hashes of consumed closed parent inputs only'))
    verify(bindings)
    return dict(root=roots['design'], design_root=roots['design'], roots=roots,
        cohorts=cohorts, controls=controls, ids=ids, labels=labels, operators=bank,
        arrays=arrays, factors=factors, parents=parents, cones=cones,
        weight_contract=receipt['source_contract'], contract=contract, closed=closed), bindings


def cohort_controls(source, c, entry):
    selected = rows(source, c); root = source['roots']['weights']
    assert entry['cohort_id'] == c['cohort_id'] and entry['records'] == len(selected)
    ap = root / entry['array_path']; jp = root / entry['record_path']
    assert sha(ap) == entry['array_sha256'] and sha(jp) == entry['record_sha256']
    rec = json.loads(jp.read_text())
    assert rec['schema'] == WEIGHT_SCHEMA and rec['source_contract'] == source['weight_contract']
    assert rec['cohort_id'] == c['cohort_id'] and rec['records'] == len(selected)
    assert rec['guide'] == c['guide'] and rec['mask'] == c['mask']
    assert rec['cohort_rows_sha256'] == c['case_rows_sha256'] and rec['ordered_case_ids_sha256'] == c['ordered_case_ids_sha256']
    assert rec['membership_occurrences'] == c['membership_occurrences']
    assert rec['original_cohort_npz_sha256'] == c['sha256'] and rec['array_file_sha256'] == sha(ap)
    assert rec['recipe'] == RECIPE and rec['residual_diagonal_prepared'] is True
    assert all(rec[k] is v for k, v in CLAIMS.items())
    with np.load(ap, allow_pickle=False) as a:
        assert a.files == ARRAYS; values = [a[k] for k in ARRAYS]
    cr, counts, weights, diagonals = values; n = len(selected)
    assert np.array_equal(cr, selected) and cr.dtype == counts.dtype == np.dtype('int64')
    assert counts.shape == (3, n) and weights.shape == diagonals.shape == (4, n)
    assert weights.dtype == diagonals.dtype == np.dtype('float64')
    assert np.all((counts >= 1) & (counts <= n))
    assert np.isfinite(weights).all() and np.isfinite(diagonals).all() and np.all(weights > 0) and np.all(diagonals > 0)
    assert np.array_equal(weights[0], np.ones(n)) and np.array_equal(diagonals[0], np.ones(n))
    assert rec['array_sha256'] == {k: array_digest(v, '<i8' if k in ARRAYS[:2] else '<f8') for k, v in zip(ARRAYS, values)}
    groups = []
    for j in range(3):
        sizes, frequencies = np.unique(counts[j], return_counts=True)
        group_count = sum((Fraction(int(f), int(s)) for s, f in zip(sizes, frequencies)), Fraction())
        assert group_count.denominator == 1; g = int(group_count); assert 0 < g <= n and n * g < 2**53
        groups.append(g)
        assert np.array_equal(weights[j + 1], n / (g * counts[j]))
        assert np.array_equal(diagonals[j + 1], (g * counts[j]) / n)
    assert rec['policy_records'] == policy_records(n, counts, weights, diagonals, groups)
    return selected, diagonals, rec


def basis(source, cohort, mode, diagonal):
    diagonal = np.asarray(diagonal)
    assert mode in MODES and diagonal.shape == (cohort['records'],)
    assert np.isfinite(diagonal).all() and np.all(diagonal > 0)
    key = cohort['cohort_id'], mode
    uniform = bool(np.array_equal(diagonal, np.ones(len(diagonal))))
    cert = source['parents'][key] if uniform else source['cones'][key]
    names = validate_certificate(cert) if uniform else cert['variance_map']['retained_names']
    return dict(names=names, exact_uniform_one=uniform, certificate_sha256=digest(cert),
        route='closed_exact_uniform_named_fold' if uniform else 'closed_positive_diagonal_cone',
        source='exact' if uniform else 'cone')


def retained_operators(source, selected, mode, names):
    assert names[0] == 'residual' and names[-1] == 'species'
    return {name: source['operators']['contrast' if name == 'family_intercept' else mode, name][selected].tocsr()
        for name in names[1:-1]}


def audit_identity(contract, design_id, mode, tree, policy):
    assert mode in MODES and policy in POLICIES
    return digest([SCHEMA, contract, design_id, mode, tree, policy])
