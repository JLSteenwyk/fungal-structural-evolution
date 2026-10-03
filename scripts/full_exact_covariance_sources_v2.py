"""General identity sources; production also requires the closed first proof."""
from pathlib import Path

from full_exact_covariance_sources import load as base_load, closed_subset, rows, local_operators, runtime_caps
from reference_measurement_union_sources import verify


def load(plan, path):
    source, bindings = base_load(plan, path)
    if plan.get('previous_completed_proof'):
        root = Path('results/phylogeny/full-exact-uniform-covariance-folds-20261003-v1')
        previous = closed_subset(plan['previous_completed_proof'], 'complete_verified_full_exact_uniform_covariance_folds',
            [root / 'cohort_certificates.jsonl', root / 'receipt.json', root / 'readback.json'], bindings)
        assert previous['cohorts'] == 4340 and previous['certificates'] == 8680
        assert previous['logical_cases'] == 75188 and previous['case_row_occurrences'] == 68220240
        source['previous_completed_proof'] = previous
    verify(bindings)
    return source, bindings
