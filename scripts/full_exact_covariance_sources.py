"""Full closed cohort/operator inputs, with explicit scoped provenance checks."""
import json
import os
from pathlib import Path

import numpy as np
from scipy import sparse

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import array_digest
from reference_measurement_union_sources import verify


def closed_subset(completion_path, status, paths, bindings):
    completion_path = Path(completion_path)
    completion = json.loads(completion_path.read_text())
    assert completion['status'] == status and completion['exact_process_journals_checked'] == 2
    assert completion['scientific_eligibility'] is False
    archive_path = Path(completion['full_hash_archive'])
    assert sha(archive_path) == completion['full_hash_archive_sha256']
    archive = json.loads(archive_path.read_text())
    assert archive['status'] == status + '_archive' and len(archive['services']) == 2
    assert len(archive['source_hashes']) == completion['bound_source_hashes']
    assert all(completion[key] == value for key, value in archive['summary'].items())
    for path in paths:
        key = str(path)
        assert key in archive['source_hashes'], key
        expected = archive['source_hashes'][key]
        assert sha(path) == expected, key
        assert key not in bindings or bindings[key] == expected
        bindings[key] = expected
    bindings[str(completion_path)] = sha(completion_path)
    bindings[str(archive_path)] = completion['full_hash_archive_sha256']
    return completion


def load(plan, path):
    bindings = dict(plan['pins']); bindings[str(path)] = sha(path)
    operator_plan = json.loads(Path(plan['operator_plan']).read_text())
    design_plan = json.loads(Path(plan['design_plan']).read_text())
    operator_root = Path(operator_plan['output']); design_root = Path(design_plan['output'])
    manifest_path = operator_root / 'operator_manifest.json'
    cohort_path = design_root / 'cohort_manifest.json'
    manifest = json.loads(manifest_path.read_text()); cohorts = json.loads(cohort_path.read_text())
    operator_paths = [Path(plan['operator_plan']), manifest_path, operator_root / 'case_ids.json', operator_root / 'block_labels.npy']
    operator_paths += [operator_root / row['path'] for row in manifest]
    operator_completion = closed_subset(plan['operator_completion'], 'complete_verified_full_entity_operator_bank', operator_paths, bindings)
    design_paths = [Path(plan['design_plan']), cohort_path] + [design_root / row['path'] for row in cohorts]
    design_completion = closed_subset(plan['design_completion'], 'complete_verified_full_expanded_model_designs', design_paths, bindings)
    assert len(cohorts) == design_completion['unique_cohorts'] == plan['expected']['cohorts'] == 4340
    assert len({row['cohort_id'] for row in cohorts}) == len(cohorts)
    assert [row['cohort_id'] for row in cohorts] == sorted(row['cohort_id'] for row in cohorts)
    ids = json.loads((operator_root / 'case_ids.json').read_text()); labels = np.load(operator_root / 'block_labels.npy', allow_pickle=False)
    assert len(ids) == len(set(ids)) == operator_completion['logical_cases'] == design_completion['logical_cases'] == 75188
    assert labels.shape == (len(ids),)
    expected = {(mode, kind) for mode in ['signed', 'unsigned'] for kind in ['target_node', 'background_node', 'model_pair', 'family', 'gene', 'model']} | {('contrast', 'family_intercept')}
    assert len(manifest) == 13 and {(row['mode'], row['kind']) for row in manifest} == expected
    operators = {}
    for row in manifest:
        value = sparse.load_npz(operator_root / row['path']).tocsr()
        assert bindings[str(operator_root / row['path'])] == row['sha256']
        assert list(value.shape) == row['shape'] and value.nnz == row['nnz'] and value.shape[0] == len(ids)
        assert value.has_canonical_format and np.isfinite(value.data).all()
        operators[row['mode'], row['kind']] = value
    _, counts = np.unique(labels, return_counts=True)
    verify(bindings)
    return dict(ids=ids, labels=labels, operators=operators, cohorts=cohorts, design_root=design_root,
        full_family_squared_rows=int(counts @ counts), operator_completion=operator_completion,
        design_completion=design_completion), bindings


def rows(source, cohort):
    path = source['design_root'] / cohort['path']
    assert sha(path) == cohort['sha256']
    with np.load(path, allow_pickle=False) as arrays:
        assert arrays.files == ['case_rows']; result = arrays['case_rows']
    assert result.dtype == np.dtype('int64') and result.shape == (cohort['records'],)
    assert len(np.unique(result)) == len(result) and len(result) > 0
    assert np.all((result >= 0) & (result < len(source['ids'])))
    assert array_digest(result, '<i8') == cohort['case_rows_sha256']
    ids = [source['ids'][index] for index in result]
    assert ids == sorted(ids) and array_digest(ids, 'S64') == cohort['ordered_case_ids_sha256']
    return result


def local_operators(source, selected, mode):
    operators = {kind: source['operators'][mode, kind][selected].tocsr()
        for kind in ['target_node', 'background_node', 'model_pair', 'gene', 'model', 'family']}
    family = source['operators']['contrast', 'family_intercept'][selected].tocsr()
    return operators, family


def runtime_caps():
    group = next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    caps = {key: (root / key).read_text().strip() for key in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert caps == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
    assert all(os.environ.get(name) == '1' for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'])
    return caps
