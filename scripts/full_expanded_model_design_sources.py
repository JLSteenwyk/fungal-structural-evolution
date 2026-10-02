"""Closed model-input I/O, schemas and exact recipe identities; no rank arithmetic."""
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
from background_measurement_union_sources import closed_source
from full_expanded_model_input_sources import ORDERS, MASKS, DISTANCE, IDENTITY, NUISANCE, STRATUM, schema
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

AXES = {'gene_distance': DISTANCE, 'aligned_identity': IDENTITY}
OUTCOMES = ['rmsd_delta', 'native_tm_dissimilarity_delta']
DEGREES = [1, 2, 3]
SETTING_FIELDS = STRATUM + ['outcome', 'sequence_axis', 'degree', 'cohort_id', 'design_id',
    'fit_input_id', 'records', 'disposition', 'nominal_tree_fits']
SUMMARY_FIELDS = ['logical_cases', 'selected_records', 'input_setting_rows', 'model_setting_rows',
    'unique_cohorts', 'unique_designs', 'unique_fit_inputs', 'nominal_tree_setting_fits',
    'unique_tree_fit_inputs', 'design_status_counts', 'fit_input_status_counts',
    'setting_status_counts', 'largest_cohort', 'cohort_member_occurrences', 'trees']


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def array_digest(a, dtype):
    original = np.asarray(a)
    a = np.ascontiguousarray(original, dtype=dtype)
    if a.dtype.kind == 'S':
        assert all(len(str(s).encode()) <= a.dtype.itemsize for s in original.reshape(-1) if not isinstance(s, bytes)), 'Identity would be truncated'
        assert all(len(s) <= a.dtype.itemsize for s in original.reshape(-1) if isinstance(s, bytes)), 'Identity would be truncated'
    if a.dtype.kind == 'f':
        a = a.copy(); a[a == 0] = 0
    return hashlib.sha256(a.tobytes()).hexdigest()


def table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def cohort_id(guide, mask, ordered_case_ids_sha256):
    return digest(['expanded-model-cohort-v1', guide, mask, ordered_case_ids_sha256])


def design_id(cohort, order, axis, degree, source_contract):
    return digest(['expanded-model-design-v1', cohort, order, axis, degree,
        ['intercept', *AXES[axis][:degree], *NUISANCE], source_contract])


def fit_id(design, outcome, response_digest):
    return digest(['expanded-model-fit-input-v1', design, outcome, response_digest])


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path)
    c = closed_source(plan['inputs_completion'], 'complete_verified_full_expanded_model_inputs',
        'complete_verified_full_expanded_model_inputs_archive', 2, bindings)
    config = json.loads(Path(plan['inputs_plan']).read_text()); bind(bindings, plan['inputs_plan'])
    root = Path(config['output']); rp = root / 'receipt.json'
    assert c['producer_receipt'] == str(rp) and bindings[str(rp)] == c['producer_receipt_sha256'] == sha(rp)
    receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_full_expanded_model_inputs_pending_independent_readback'
    assert receipt['plan_sha256'] == sha(plan['inputs_plan']) and receipt['scientific_eligibility'] is False
    for key in ['logical_cases', 'selected_records']:
        assert c[key] == receipt[key] == plan['expected'][key]
    assert c['settings'] == plan['expected']['input_setting_rows']
    assert c['future_model_setting_rows'] == plan['expected']['model_setting_rows']
    for name, d in receipt['artifacts'].items(): bind(bindings, root / name, d)
    del receipt
    cp = Path(config['cases_plan']); bind(bindings, cp)
    case_config = json.loads(cp.read_text()); case_root = Path(case_config['output'])
    matching = json.loads(Path(case_config['matching_plan']).read_text())
    scenarios_path = Path(matching['selection']) / 'scenarios.json'; assert str(scenarios_path) in bindings
    scenarios = [r['scenario_id'] for r in json.loads(scenarios_path.read_text())]
    assert len(scenarios) == len(set(scenarios)) == config['expected']['scenarios']
    cases = table(case_root / 'case_index.tsv.gz')
    assert len(cases) == len({r['case_id'] for r in cases}) == c['logical_cases']
    vp = Path(config['covariance_plan']); bind(bindings, vp)
    covariance = json.loads(vp.read_text()); covroot = Path(covariance['output'])
    cov = {r['case_id']: r for r in table(covroot / 'case_covariance_index.tsv.gz')}
    assert set(cov) == {r['case_id'] for r in cases}
    vc = json.loads(Path(config['covariance_completion']).read_text())
    assert vc['status'] == 'complete_verified_full_expanded_covariance' and vc['scientific_eligibility'] is False
    assert vc['trees'] == covariance['trees'] == plan['trees']
    parts = json.loads((root / 'partition_manifest.json').read_text())
    assert [(p['mask'], p['order_contrast']) for p in parts] == [(m, o) for m in MASKS for o in ORDERS]
    arrays = {}
    names = ['input_id', 'case_id', 'case_row', 'numerical_usable', 'joint_mask_pass_bits', 'joint_both_pass_bits',
        *OUTCOMES, *DISTANCE, *IDENTITY, *NUISANCE]
    for p in parts:
        fp = root / p['path']; bind(bindings, fp, p['sha256']); pf = pq.ParquetFile(fp)
        assert pf.schema_arrow == schema() and pf.metadata.num_rows == p['rows'] == len(cases)
        t = pf.read(columns=names)
        assert t['case_id'].to_pylist() == [r['case_id'] for r in cases]
        assert t['case_row'].to_pylist() == list(range(len(cases)))
        assert all(len(s.encode()) <= 64 for k in ['input_id','case_id'] for s in t[k].to_pylist())
        arrays[p['mask'], p['order_contrast']] = {
            k: np.asarray(t[k].to_pylist(), dtype='S64' if k in ['input_id', 'case_id'] else np.int64 if k in names[2:6] else np.float64)
            for k in names}
    counts = table(root / 'setting_counts.tsv')
    assert len(counts) == c['settings'] and len({tuple(r[k] for k in STRATUM) for r in counts}) == len(counts)
    assert {tuple(r[k] for k in STRATUM) for r in counts} == set(itertools.product(
        config['guides'], config['policies'], scenarios, MASKS, ['mask', 'both_masks'],
        [s['id'] for s in config['screens']], ORDERS))
    links = case_root / 'selection_case_links.tsv.gz'; assert str(links) in bindings
    for p in [case_root / 'case_index.tsv.gz', covroot / 'case_covariance_index.tsv.gz', config['covariance_completion']]:
        assert str(p) in bindings
    contract = digest(dict(input_plan=sha(plan['inputs_plan']), input_completion=sha(plan['inputs_completion']),
        input_producer=c['producer_receipt_sha256'], covariance_completion=sha(config['covariance_completion']),
        covariance_plan=sha(vp), schema='expanded-model-design-v1'))
    verify(bindings)
    return dict(cases=cases, cov=cov, config=config, root=root, arrays=arrays, counts=counts,
        selections=links, scenarios=scenarios, contract=contract), bindings
