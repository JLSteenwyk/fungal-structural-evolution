"""Full closed structural fits and all-order summaries for contrast sensitivity."""
import json
from pathlib import Path
from full_triad_context_geometry_sources import load_sources as closed_groups
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha

MASKS = ['full', 'plddt70', 'both_masks']
CORES = ['reference_common', 'cycle_consistent', 'both_cores']
DIRECTIONS = ['positive', 'negative', 'within_numerical_tolerance', 'sign_uncertain',
              'unavailable', 'nonunique_fit']
QUALIFIED = DIRECTIONS + ['excluded_by_quality']
SUMMARY_FIELDS = ['measured_triads', 'source_fit_rows', 'source_order_groups',
                  'sensitivity_groups', 'screen_decisions', 'summary_rows',
                  'direction_counts', 'qualified_direction_counts', 'maximum_contrast_span']


def load_sources(plan, path):
    _, groups, _, bindings = closed_groups(plan, path)
    geometry = json.loads(Path(plan['geometry_completion']).read_text())
    assert geometry['status'] == 'complete_verified_full_triad_same_residue_geometry'
    assert geometry['exact_process_journals_checked'] == 2
    rp = Path(geometry['producer_receipt'])
    assert str(rp) in bindings and bindings[str(rp)] == geometry['producer_receipt_sha256']
    raw = rp.parent / 'common_residue_fits.tsv.gz'
    assert str(raw) in bindings
    work = json.loads(Path(plan['triad_work_plan']).read_text())
    catalog = Path(work['output']) / 'ordered_model_triads.jsonl'
    assert str(catalog) in bindings
    with catalog.open() as f:
        triads = [r for line in f if (r := json.loads(line))['source_design_ready_links']]
    assert len(triads) == plan['expected']['measured_triads'] == 27056
    assert geometry['fit_rows'] == plan['expected']['source_fit_rows'] == 865792
    assert plan['expected']['source_order_groups'] == 108224
    assert plan['expected']['sensitivity_groups'] == 243504
    assert plan['expected']['screen_decisions'] == 1461024
    assert plan['expected']['summary_rows'] == 378
    assert plan['contrast_numerical_tolerance_angstrom'] == 1e-9
    bind(bindings, plan['geometry_completion'])
    verify(bindings)
    return raw, groups, triads, bindings


def selected_keys(mask, core):
    return [(m, c) for m in MASKS[:2] for c in CORES[:2]
            if (mask == 'both_masks' or mask == m) and (core == 'both_cores' or core == c)]
