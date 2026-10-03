"""Fresh hashes for every consumed original case/cohort, with scoped closures."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset, rows
from full_expanded_model_design_sources import cohort_id, digest
from inverse_reuse_weight_controls import SCHEMA
from reference_measurement_union_sources import bind, verify

SOURCES = {'cases': 'complete_verified_full_matching_logical_case_index',
           'covariance': 'complete_verified_full_expanded_covariance',
           'operator': 'complete_verified_full_entity_operator_bank',
           'design': 'complete_verified_full_expanded_model_designs',
           'cone': 'complete_verified_full_positive_diagonal_covariance_cones_v1'}
SHARED = ['case_id', 'physical_case_id', 'target_id', 'background_id', 'guide',
          'target_family', 'background_family', 'selection_records',
          'target_same_model', 'background_same_model']


def read_columns(path, names):
    with gzip.open(path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t'); assert set(names) <= set(reader.fieldnames)
        return [{k: r[k] for k in names} for r in reader]


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path); verify(bindings)
    configs = {k: json.loads(Path(plan[k + '_plan']).read_text()) for k in SOURCES}
    roots = {k: Path(v['output']) for k, v in configs.items()}
    assert configs['operator']['covariance_plan'] == plan['covariance_plan']
    assert configs['covariance']['case_plan'] == plan['cases_plan']
    assert configs['covariance']['case_completion'] == plan['cases_completion']
    assert configs['operator']['covariance_completion'] == plan['covariance_completion']
    input_path = Path(configs['design']['inputs_plan'])
    inputs = json.loads(input_path.read_text())
    assert inputs['cases_plan'] == plan['cases_plan']
    assert inputs['covariance_plan'] == plan['covariance_plan']
    manifest = roots['design'] / 'cohort_manifest.json'
    cohorts = json.loads(manifest.read_text())
    assert len(cohorts) == plan['expected']['cohorts']
    assert [c['cohort_id'] for c in cohorts] == sorted({c['cohort_id'] for c in cohorts})
    paths = {
        'cases': [roots['cases'] / 'case_index.tsv.gz'],
        'covariance': [roots['covariance'] / 'case_covariance_index.tsv.gz'],
        'operator': [roots['operator'] / 'case_ids.json'],
        'design': [manifest, input_path, *[roots['design'] / c['path'] for c in cohorts]],
        'cone': [roots['cone'] / 'cohort_cones.jsonl'],
    }
    closed = {}
    for key in SOURCES:
        closed[key] = closed_subset(plan[key + '_completion'], SOURCES[key],
                                   [Path(plan[key + '_plan']), *paths[key]], bindings)
    ids = json.loads(paths['operator'][0].read_text()); n = plan['expected']['logical_cases']
    assert len(ids) == len(set(ids)) == n
    assert all(closed[k]['logical_cases'] == n for k in SOURCES)
    assert closed['design']['unique_cohorts'] == closed['cone']['cohorts'] == len(cohorts)
    cases = read_columns(paths['cases'][0], SHARED + ['background_pair_key'])
    cov = read_columns(paths['covariance'][0], SHARED + ['family_component'])
    assert [r['case_id'] for r in cases] == [r['case_id'] for r in cov] == ids
    columns = {'background_node': [], 'background_pair': [], 'family_component': []}
    guides = []
    for a, b in zip(cases, cov):
        assert all(a[k] == b[k] for k in SHARED)
        for value in [a['background_id'], a['background_pair_key'], b['family_component']]:
            assert len(value) == 64 and all(ch in '0123456789abcdef' for ch in value)
        columns['background_node'].append(a['background_id'])
        columns['background_pair'].append(a['background_pair_key'])
        columns['family_component'].append(b['family_component'])
        guides.append(a['guide'])
    del cases, cov
    source = dict(ids=ids, design_root=roots['design'], cohorts=cohorts,
                  keys={k: np.asarray(v, dtype='S64') for k, v in columns.items()}, guides=guides)
    by_id = {c['cohort_id']: c for c in cohorts}; seen = set(); cone_rows = 0
    with paths['cone'][0].open() as f:
        for line in f:
            row = json.loads(line); c = by_id[row['cohort_id']]
            mode = row['loading_mode']; assert mode in ['signed', 'unsigned']
            key = c['cohort_id'], mode; assert key not in seen; seen.add(key)
            assert row['records'] == c['records']
            assert row['cohort_rows_sha256'] == c['case_rows_sha256']
            assert row['ordered_case_ids_sha256'] == c['ordered_case_ids_sha256']
            assert row['scientific_eligibility'] is False
            assert row['variance_map']['residual_diagonal_prepared'] is False
            assert row['variance_map']['raw_reml_basis_qualification_complete'] is False
            cone_rows += row['records']
    assert seen == {(c['cohort_id'], mode) for c in cohorts for mode in ['signed', 'unsigned']}
    assert len(seen) == closed['cone']['certificates'] == 2 * len(cohorts)
    assert cone_rows == closed['cone']['case_row_occurrences']
    occurrences = 0
    for c in cohorts:
        assert c['cohort_id'] == cohort_id(c['guide'], c['mask'], c['ordered_case_ids_sha256'])
        selected = rows(source, c)
        assert all(guides[i] == c['guide'] for i in selected)
        assert bindings[str(roots['design'] / c['path'])] == c['sha256']
        occurrences += len(selected)
    assert occurrences == closed['design']['cohort_member_occurrences'] == plan['expected']['case_row_occurrences']
    assert cone_rows == 2 * occurrences
    verify(bindings)
    contract = digest(dict(schema=SCHEMA, source_hashes=bindings,
                          group_definition='cohort-local logical cases; original background node, versioned pair, connected family component'))
    return source, contract, bindings
