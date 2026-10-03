#!/usr/bin/env python3
"""Certify named dependencies over every cohort and both loading modes."""
import argparse
from collections import Counter
import fcntl
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from covariance_exact_folds_v2 import certificate, variance_map
from full_exact_covariance_sources_v2 import load, rows, local_operators, runtime_caps
from reference_measurement_union_sources import verify

STATUS = 'complete_full_exact_uniform_covariance_fold_certificates_pending_readback_v2'
SUMMARY_FIELDS = ['logical_cases', 'cohorts', 'loading_modes', 'certificates', 'case_row_occurrences',
                  'retained_basis_counts', 'exact_relation_counts', 'failed_relation_counts']


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text()); assert plan.get('previous_completed_proof')
    caps = runtime_caps()
    source, bindings = load(plan, args.plan); root = Path(plan['output'])
    root.mkdir(exist_ok=False); lock = (root / 'stage.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    counts = Counter(); positive = Counter(); negative = Counter(); total = occurrences = 0
    with (root / 'cohort_certificates.jsonl').open('x') as handle:
        for index, cohort in enumerate(source['cohorts']):
            selected = rows(source, cohort)
            for mode in ['signed', 'unsigned']:
                operators, family = local_operators(source, selected, mode)
                proof = certificate(operators, family, mode, selected)
                mapping = variance_map(proof['relations'], mode)
                record = dict(cohort_id=cohort['cohort_id'], cohort_rows_sha256=cohort['case_rows_sha256'],
                    ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'], certificate=proof, variance_map=mapping,
                    raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
                handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + '\n')
                total += 1; occurrences += len(selected); counts[str(len(mapping['retained_names']))] += 1
                for name, result in proof['relations'].items():
                    (positive if result['exact'] else negative)[mode + ':' + name] += 1
            if (index + 1) % 20 == 0: print('exact_covariance_cohorts', index + 1, '/4340', flush=True)
    assert total == 8680
    summary = dict(logical_cases=75188, cohorts=4340, loading_modes=['signed', 'unsigned'], certificates=total,
        case_row_occurrences=occurrences, retained_basis_counts=dict(counts),
        exact_relation_counts=dict(positive), failed_relation_counts=dict(negative))
    verify(bindings)
    result = dict(status=STATUS, plan_sha256=sha(args.plan), **summary, source_hashes=bindings,
        actual_cgroup_limits=caps,
        artifacts={'cohort_certificates.jsonl': sha(root / 'cohort_certificates.jsonl')},
        scientific_eligibility=False, scope=plan['scope'])
    with (root / 'receipt.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__': main()
