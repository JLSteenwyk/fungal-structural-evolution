"""Bind the dedicated closed graph/covariate proof and its matching lineage."""
import json
from pathlib import Path
from reference_measurement_union_sources import bind
from run_ortholog_pair_guide_comparison import sha


def load_graph(plan, bindings, matching_receipt):
    completion = Path(plan['graph_completion'])
    closed = json.loads(completion.read_text()); bind(bindings, completion)
    assert closed['status'] == 'complete_verified_expanded_background_graph_covariates'
    assert len(closed['services']) == 4 and closed['scientific_eligibility'] is False
    for path, digest in closed['source_hashes'].items(): bind(bindings, path, digest)
    graph = Path(plan['graph']); covariates = Path(plan['covariates'])
    gp, cp = graph / 'receipt.json', covariates / 'receipt.json'
    gap, cap = Path(plan['graph_readback']), Path(plan['covariate_readback'])
    assert all(str(p) in closed['source_hashes'] for p in [gp, cp, gap, cap])
    gr, cr, ga, ca = [json.loads(p.read_text()) for p in [gp, cp, gap, cap]]
    assert gr['status'] == 'complete_background_match_graph_pending_readback'
    assert cr['status'] == 'complete_background_match_covariates_pending_readback'
    assert ga['status'] == 'passed_full_background_match_graph_readback'
    assert ca['status'] == 'passed_full_background_match_covariate_readback'
    assert ga['producer_receipt_sha256'] == cr['graph_receipt_sha256'] == sha(gp)
    assert cr['graph_readback_sha256'] == sha(gap)
    assert ca['producer_receipt_sha256'] == matching_receipt['source_covariate_receipt_sha256'] == sha(cp)
    summary = closed['summary']
    assert gr['target_nodes'] == ga['target_nodes'] == summary['duplicate_target_links']
    assert gr['background_nodes'] == ga['background_nodes'] == summary['background_nodes']
    assert gr['target_policy_dispositions'] == ga['target_policy_dispositions'] == summary['target_policy_dispositions']
    assert gr['edges'] == ga['edges'] == cr['edges'] == ca['edges'] == matching_receipt['source_edges'] == summary['eligible_edges']
    assert ga['support_rows_checked'] == summary['architecture_support_rows_checked']
    for name, digest in gr['artifacts'].items():
        assert closed['source_hashes'].get(str(graph / name)) == digest
        bind(bindings, graph / name, digest)
    return gr
