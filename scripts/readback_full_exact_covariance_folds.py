#!/usr/bin/env python3
"""Rebuild every exact certificate from independent CSC column outer sums."""
import argparse
from collections import Counter
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from covariance_exact_folds import independent_certificate, independent_variance_map
from full_exact_covariance_sources import load, rows, local_operators, runtime_caps
from reference_measurement_union_sources import verify
from run_full_exact_covariance_folds import STATUS as PRODUCER_STATUS, SUMMARY_FIELDS

STATUS = 'passed_full_exact_uniform_covariance_fold_independent_readback'


def check_record(record, cohort, mode, operators, family, selected):
    expected = independent_certificate(operators, family, mode, selected)
    mapping = independent_variance_map(expected['relations'], mode)
    reconstructed = dict(cohort_id=cohort['cohort_id'], cohort_rows_sha256=cohort['case_rows_sha256'],
        ordered_case_ids_sha256=cohort['ordered_case_ids_sha256'], certificate=expected, variance_map=mapping,
        raw_reml_basis_qualification_complete=False, scientific_eligibility=False)
    assert record == reconstructed
    return expected, mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args(); plan = json.loads(args.plan.read_text()); source, bindings = load(plan, args.plan)
    caps = runtime_caps()
    root = Path(plan['output']); receipt_path = root / 'receipt.json'; producer = json.loads(receipt_path.read_text())
    assert producer['status'] == PRODUCER_STATUS and producer['plan_sha256'] == sha(args.plan)
    assert not (root / 'readback.json').exists()
    verify(producer['source_hashes'])
    assert sha(root / 'cohort_certificates.jsonl') == producer['artifacts']['cohort_certificates.jsonl']
    counts = Counter(); positive = Counter(); negative = Counter(); total = occurrences = 0
    with (root / 'cohort_certificates.jsonl').open() as handle:
        for index, cohort in enumerate(source['cohorts']):
            selected = rows(source, cohort)
            for mode in ['signed', 'unsigned']:
                record = json.loads(next(handle)); operators, family = local_operators(source, selected, mode)
                proof, mapping = check_record(record, cohort, mode, operators, family, selected)
                total += 1; occurrences += len(selected); counts[str(len(mapping['retained_names']))] += 1
                for name, result in proof['relations'].items():
                    (positive if result['exact'] else negative)[mode + ':' + name] += 1
            if (index + 1) % 20 == 0: print('independent_exact_covariance_cohorts', index + 1, '/4340', flush=True)
        assert next(handle, None) is None, 'Foreign certificates retained outside full original grid'
    assert total == 8680
    summary = dict(logical_cases=75188, cohorts=4340, loading_modes=['signed', 'unsigned'], certificates=total,
        case_row_occurrences=occurrences, retained_basis_counts=dict(counts),
        exact_relation_counts=dict(positive), failed_relation_counts=dict(negative))
    assert all(summary[key] == producer[key] for key in SUMMARY_FIELDS)
    verify(bindings); verify(producer['source_hashes'])
    assert sha(root / 'cohort_certificates.jsonl') == producer['artifacts']['cohort_certificates.jsonl']
    bindings[str(receipt_path)] = sha(receipt_path)
    result = dict(status=STATUS, plan_sha256=sha(args.plan), producer_receipt_sha256=sha(receipt_path),
        actual_cgroup_limits=caps,
        **summary, source_hashes=bindings, scientific_eligibility=False, scope=plan['scope'])
    with (root / 'readback.json').open('x') as handle: handle.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__': main()
