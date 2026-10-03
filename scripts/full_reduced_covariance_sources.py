"""Scoped closed provenance for the complete retained-kernel qualification."""
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_exact_covariance_sources import closed_subset
from full_expanded_model_design_sources import digest, SETTING_FIELDS
from full_covariance_qualification_sources import jsonl, LINK_EXTRA
from reduced_covariance_basis import SCHEMA, validate_certificate
from reference_measurement_union_sources import verify

NEW_LINK_EXTRA = LINK_EXTRA + ['source_audit_id', 'source_covariance_disposition']
SUMMARY_FIELDS = ['logical_cases', 'model_setting_rows', 'unique_cohorts', 'unique_designs',
    'audit_rows', 'setting_audit_links', 'audit_status_counts', 'link_status_counts', 'trees',
    'loading_modes', 'retained_basis_audit_counts']


def load(plan, path):
    bindings = dict(plan['pins']); bindings[str(path)] = sha(path)
    old_plan_path = Path(plan['qualification_plan']); old_plan = json.loads(old_plan_path.read_text())
    old_root = Path(old_plan['output'])
    old_paths = [old_plan_path] + [old_root / name for name in ['receipt.json', 'readback.json',
        'stage_plan.json', 'design_covariance_audits.jsonl.gz', 'setting_audit_links.tsv.gz']]
    old = closed_subset(plan['qualification_completion'], 'complete_verified_full_uniform_covariance_qualification', old_paths, bindings)
    exact_plan_path = Path(plan['exact_plan']); exact_plan = json.loads(exact_plan_path.read_text())
    exact_root = Path(exact_plan['output'])
    exact_paths = [exact_plan_path] + [exact_root / name for name in ['receipt.json', 'readback.json', 'cohort_certificates.jsonl']]
    exact = closed_subset(plan['exact_completion'], 'complete_verified_full_exact_uniform_covariance_folds_v2', exact_paths, bindings)
    expected = plan['expected']
    for field in ['logical_cases', 'model_setting_rows', 'unique_cohorts', 'unique_designs', 'audit_rows', 'setting_audit_links']:
        assert old[field] == expected[field], field
    assert exact['logical_cases'] == old['logical_cases'] and exact['cohorts'] == old['unique_cohorts']
    assert old['trees'] == plan['trees'] and old['loading_modes'] == ['signed', 'unsigned']
    assert old['audit_rows'] == old['unique_designs'] * len(old['trees']) * len(old['loading_modes'])
    assert old['setting_audit_links'] == old['model_setting_rows'] * len(old['trees']) * len(old['loading_modes'])
    certificates = {}
    for row in jsonl(exact_root / 'cohort_certificates.jsonl'):
        validate_certificate(row)
        key = (row['cohort_id'], row['certificate']['loading_mode'])
        assert key not in certificates; certificates[key] = row
    assert len(certificates) == exact['certificates'] == 2 * old['unique_cohorts']
    cohort_ids = sorted({k[0] for k in certificates})
    assert set(certificates) == {(c, m) for c in cohort_ids for m in old['loading_modes']}
    contract = digest(dict(schema=SCHEMA, original_qualification_completion=sha(plan['qualification_completion']),
        exact_completion=sha(plan['exact_completion']), original_audits_sha256=bindings[str(old_root / 'design_covariance_audits.jsonl.gz')],
        certificates_sha256=bindings[str(exact_root / 'cohort_certificates.jsonl')],
        numerical_envelopes='copy_original_principal_submatrices_unchanged', residual='uniform_one'))
    verify(bindings)
    return dict(old=old, old_root=old_root, certificates=certificates, cohort_ids=cohort_ids,
        contract=contract, settings_fields=SETTING_FIELDS), bindings
