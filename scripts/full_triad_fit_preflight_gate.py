"""Require the closed full-work preflight before production geometric fitting."""
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def require_preflight(plan, bindings):
    path = plan['preflight_completion']; bind(bindings, path)
    c = json.loads(Path(path).read_text())
    assert c['status'] == 'complete_verified_full_triad_geometric_fit_work_preflight'
    assert c['exact_process_journals_checked'] == 1
    bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256'])
    archive = json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status'] == c['status'] and len(archive['services']) == 1
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    assert archive['summary']['mapping_states'] == c['mapping_states'] == plan['expected']['mask_order_states']
    assert archive['summary']['fit_rows'] == c['fit_rows'] == 2 * c['mapping_states']
    assert archive['source_hashes'][plan['mapping_completion']] == sha(plan['mapping_completion'])
    assert archive['source_hashes'][c['preflight_receipt']] == c['preflight_receipt_sha256']
    r = json.loads(Path(c['preflight_receipt']).read_text())
    assert r['status'] == 'passed_full_triad_geometric_fit_work_preflight'
    assert all(r[key] == value for key, value in archive['summary'].items())
    for source, digest in archive['source_hashes'].items(): bind(bindings, source, digest)
    verify(bindings)
