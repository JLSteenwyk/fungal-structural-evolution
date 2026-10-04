#!/usr/bin/env python3
"""Independently retain a completed native V6 role with explicit special values."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from baliphy_joint_sampler_scalar_v6 import inspect, REVIEW
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--chain', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    plan = json.loads(a.plan.read_text()); verify(plan['pins'])
    jobs = json.loads(Path(plan['jobs']).read_text())
    selected = [j for j in jobs if j['chain']['chain_id'] == a.chain]; assert len(selected) == 1
    job = selected[0]; verify(job['config']['pins'])
    root = Path(plan['output']); row_path = root/'chains'/(a.chain+'.json')
    row = json.loads(row_path.read_text()); assert row['status'] == REVIEW
    native_path = Path(row['native_receipt']); native = json.loads(native_path.read_text())
    assert row['native_receipt_sha256'] == sha(native_path)
    assert native['exit_code'] == row['exit_code'] == 0
    rebuilt = inspect(job, native_path, sha(a.plan), plan['mapping'], root/'frames'/a.chain,
        allow_export_creation=False)
    assert rebuilt == row and not rebuilt['scalar_integrity_accepted'] and not rebuilt['joint_frames']
    assert not (root/'frames'/a.chain).exists()
    tsv = Path(row['sample_audit']['scalar_log']); scalar_json = tsv.with_name('C1.log.json')
    records = [json.loads(line) for line in scalar_json.read_text().splitlines()][1:]
    table = list(csv.DictReader(tsv.read_text().splitlines(), delimiter='\t'))
    assert len(records) == len(table) == 21
    observations = []
    for review in row['scalar_v6_audit']['nonfinite_reviews']:
        iteration = review['iteration']; value = records[iteration][{'parameters': 'parameters//', 'context': 'statistics//'}[review['section']]]
        for key in review['path']: value = value[key]
        assert value == '__project_scalar_v6__:' + review['kind']
        observations.append(dict(**review, encoded_json_value=value,
            tsv_parameter_token=table[iteration].get(review['path'][-1])))
    assert observations or row['scalar_v6_audit']['literal_null_iterations']
    pins = dict(plan['pins']); modules = project_sources(pins, [Path(__file__)])
    for path in [a.plan, row_path, native_path, scalar_json, tsv]: bind(pins, path)
    for name, digest in native['artifacts'].items(): bind(pins, native_path.parent/name, digest)
    verify(pins)
    result = dict(status='verified_completed_native_scalar_v6_special_value_role_retained',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan=str(a.plan), plan_sha256=sha(a.plan),
        chain_id=a.chain, family=row['family'], prior_label=row['prior_label'], seed=row['seed'],
        native_exit_code=0, independent_role_reconstruction_matched=True,
        scalar_rows_checked=21, special_value_observations=observations,
        literal_null_iterations=row['scalar_v6_audit']['literal_null_iterations'],
        admitted_ancestral_arrays=0, transitive_project_source_modules=len(modules), source_hashes=pins,
        posterior_qualified=False, scientific_eligibility=False, gpu=False,
        scope='One completed actual original V6 role, independent strict scalar audit and exact whole '
              'saved disposition reconstruction without creating exports. Explicit special-value native '
              'JSON tokens retained alongside TSV tokens; no claim about the proposal/model cause. '
              'Unresolved role stays in complete grid accounting. Not a new native attempt, retry, '
              'prior adjustment, full-grid readback, posterior or biological acceptance.')
    with a.output.open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__': main()
