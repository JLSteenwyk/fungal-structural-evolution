"""Closed complete exact operator certificates for diagonal-independent cones."""
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from background_measurement_union_sources import closed_source
from full_expanded_model_design_sources import digest
from nonuniform_covariance_cone import SCHEMA, validate_parent
from reference_measurement_union_sources import bind, verify


def load(plan, path):
    bindings = dict(plan['pins']); bind(bindings, path)
    status = 'complete_verified_full_exact_uniform_covariance_folds_v2'
    completion = closed_source(plan['parent_completion'], status, status + '_archive', 2, bindings)
    assert completion['scientific_eligibility'] is False
    config = json.loads(Path(plan['parent_plan']).read_text()); bind(bindings, plan['parent_plan'])
    cert_path = Path(config['output']) / 'cohort_certificates.jsonl'
    assert str(cert_path) in bindings
    verify(bindings)
    with cert_path.open() as handle:
        parents = [json.loads(line) for line in handle]
    assert len(parents) == completion['certificates'] == plan['expected']['certificates']
    assert completion['logical_cases'] == plan['expected']['logical_cases']
    assert completion['cohorts'] == plan['expected']['cohorts']
    seen = {}; occurrences = 0
    for parent in parents:
        c = validate_parent(parent); key = parent['cohort_id'], c['loading_mode']
        assert key not in seen; seen[key] = parent
        occurrences += c['records']
    cohort_ids = sorted({key[0] for key in seen})
    assert len(cohort_ids) == completion['cohorts']
    assert set(seen) == {(cid, mode) for cid in cohort_ids for mode in ['signed', 'unsigned']}
    for cid in cohort_ids:
        a, b = [seen[cid, mode] for mode in ['signed', 'unsigned']]
        for field in ['cohort_rows_sha256', 'ordered_case_ids_sha256']:
            assert a[field] == b[field]
        assert a['certificate']['records'] == b['certificate']['records']
    if 'case_row_occurrences' in completion:
        assert occurrences == completion['case_row_occurrences']
    contract = digest(dict(schema=SCHEMA, parent_completion=sha(plan['parent_completion']),
        parent_plan=sha(plan['parent_plan']), certificates=sha(cert_path), pins=plan['pins']))
    return parents, completion, contract, bindings
