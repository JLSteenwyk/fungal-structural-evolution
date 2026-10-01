"""Closed full identical-residue fits and complete logical work-design source I/O."""
import json
from pathlib import Path
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

DEFINITIONS = ['reference_common', 'cycle_consistent']
NUMERIC_FIELDS = ['common_residues', 'coverage_a', 'coverage_b', 'coverage_reference',
                  'rmsd_ab', 'rmsd_ar', 'rmsd_br', 'rmsd_ar_minus_br',
                  'sequence_identity_ab', 'sequence_identity_ar', 'sequence_identity_br',
                  'mean_plddt_a', 'mean_plddt_b', 'mean_plddt_reference',
                  'joint_plddt70_fraction', 'relative_curvature_ab',
                  'relative_curvature_ar', 'relative_curvature_br']
SUMMARY_FIELDS = ['correspondence_work_triads', 'fit_rows', 'robustness_groups',
                  'counts', 'all_order_screen_pass_counts', 'any_order_screen_pass_counts',
                  'maximum_metric_spans']


def load_sources(plan, plan_path):
    bindings = dict(plan['pins']); bind(bindings, plan_path)
    cp = plan['geometry_completion']; c = json.loads(Path(cp).read_text())
    assert c['status'] == 'complete_verified_full_triad_same_residue_geometry' and c['exact_process_journals_checked'] == 2
    bind(bindings, cp); bind(bindings, c['full_hash_archive'], c['full_hash_archive_sha256'])
    archive = json.loads(Path(c['full_hash_archive']).read_text())
    assert archive['status'] == 'complete_verified_full_triad_same_residue_fit_archive' and len(archive['services']) == 2
    assert len(archive['source_hashes']) == c['bound_source_hashes']
    for path, digest in archive['source_hashes'].items(): bind(bindings, path, digest)
    rp = Path(c['producer_receipt']); ap = Path(c['independent_readback'])
    bind(bindings, rp, c['producer_receipt_sha256']); bind(bindings, ap, c['independent_readback_sha256'])
    r, a = [json.loads(p.read_text()) for p in [rp, ap]]
    assert r['status'] == 'complete_full_triad_same_residue_fits_pending_independent_readback'
    assert a['status'] == 'passed_full_triad_same_residue_quaternion_readback' and a['producer_receipt_sha256'] == sha(rp)
    assert r['plan_sha256'] == a['plan_sha256'] == sha(plan['geometry_plan'])
    for key in ['correspondence_work_triads', 'mapping_states', 'fit_rows']:
        assert r[key] == a[key] == c[key] == plan['expected'][key]
    config = json.loads(Path(plan['geometry_plan']).read_text()); assert config['screens'] == plan['screens']
    for name, digest in r['artifacts'].items(): bind(bindings, rp.parent / name, digest)
    tp = plan['triad_work_plan']; tc = json.loads(Path(plan['triad_work_completion']).read_text())
    assert tc['status'] == 'complete_verified_full_reference_triad_work_design' and len(tc['services']) == 2
    for path, digest in tc['source_hashes'].items(): bind(bindings, path, digest)
    tr = Path(json.loads(Path(tp).read_text())['output']); catalog = []
    with (tr / 'ordered_model_triads.jsonl').open() as f:
        for line in f: catalog.append(json.loads(line))
    assert len(catalog) == tc['summary']['unique_ordered_model_triads'] == 31235
    triads = [t for t in catalog if t['source_design_ready_links']]
    assert len(triads) == plan['expected']['correspondence_work_triads']
    assert all(t['distinct_versioned_models'] == t['distinct_model_ids'] == 3 for t in triads)
    bind(bindings, plan['geometry_plan']); bind(bindings, tp); bind(bindings, plan['triad_work_completion'])
    verify(bindings)
    return rp.parent / 'common_residue_fits.tsv.gz', triads, tr / 'context_triad_design.jsonl.gz', tc, bindings
