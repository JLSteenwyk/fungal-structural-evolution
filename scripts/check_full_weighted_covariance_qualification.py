#!/usr/bin/env python3
"""Complete synthetic four-control numerical/export/readback and corruption grid."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from check_full_inverse_reuse_weights import closure, write
from check_full_weighted_covariance_sources_v4 import setup, rejected
from full_weighted_covariance_qualification import run, SUMMARY
from reference_measurement_union_sources import verify
from run_full_weighted_covariance_source_census_v3 import run as census_run, SUMMARY as CENSUS_SUMMARY


def fixture(root):
    pp = setup(root); source_plan = json.loads(pp.read_text())
    producer = census_run(pp); reader = census_run(pp, True); sr = Path(source_plan['output'])
    cp = root / 'census.completed.json'
    paths = set(map(Path, reader['source_hashes'])) | {pp, sr / 'receipt.json', sr / 'readback.json', sr / 'cohort_source_census.jsonl.gz'}
    closure(cp, 'complete_verified_full_four_control_covariance_source_census_v1',
        sorted(paths), {k:reader[k] for k in CENSUS_SUMMARY})
    plan = dict(source_plan); plan['output'] = str(root / 'numerical')
    plan.update(source_census_plan=str(pp), source_census_completion=str(cp),
        resources=dict(cpus=2, memory_gib=16),
        scope='Complete declared synthetic five original nonempty cohorts/150designs/600settings, all four controls/two modes/five trees:6000audit records and24000original setting links. Parent source/journal closures synthetic. No real numerical qualification, fit or biological pilot.')
    modules = ['full_weighted_covariance_qualification', 'check_full_weighted_covariance_qualification',
        'independent_positive_diagonal_basis_context', 'independent_positive_diagonal_kernel_products',
        'positive_diagonal_basis_context', 'covariance_basis_context', 'covariance_basis_audit',
        'readback_full_covariance_qualification', 'reduced_covariance_basis', 'full_weighted_covariance_sources_v2',
        'check_full_weighted_covariance_sources_v4', 'run_full_weighted_covariance_source_census_v3',
        'check_full_covariance_qualification', 'check_full_entity_operators', 'prepare_full_entity_operators',
        'readback_full_entity_operators', 'full_entity_operator_sources', 'full_exact_covariance_sources',
        'full_expanded_model_design_sources', 'full_expanded_model_input_sources', 'nonuniform_covariance_cone',
        'covariance_exact_folds_v2', 'inverse_reuse_weight_controls', 'run_full_inverse_reuse_weights',
        'reference_measurement_union_sources', 'run_positive_diagonal_kernel_software_stage']
    plan['pins'] = {'scripts/' + m + '.py':sha('scripts/' + m + '.py') for m in modules}
    np = root / 'numerical.plan.json'; write(np, plan); return np


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True); a = p.parse_args()
    a.output.mkdir(exist_ok=False); assert not a.receipt.exists(); pp = fixture(a.output / 'fixture')
    plan = json.loads(pp.read_text()); root = Path(plan['output'])
    producer = run(pp); reader = run(pp, True)
    assert (producer['cohorts'], producer['designs'], producer['settings'], producer['numerical_audit_rows'], producer['setting_audit_links']) == (5, 150, 600, 6000, 24000)
    assert all(reader[k] == producer[k] for k in SUMMARY)
    assert reader['producer_receipt_sha256'] == sha(root / 'receipt.json')
    assert producer['working_model_fits_computed'] == 0 and producer['nonuniform_weighting_accepted'] is False
    rejected(lambda: run(pp)); rejected(lambda: run(pp, True))
    rp = root / 'receipt.json'; rb = root / 'readback.json'; mp = root / 'cohort_manifest.json'
    manifest = json.loads(mp.read_text())
    entry = next(e for e in manifest if any(json.loads(l)['numerical_audit'] is not None
        for l in gzip.decompress((root / e['audit_file']).read_bytes()).decode().splitlines()))
    ap = root / entry['audit_file']; lp = root / entry['link_file']
    originals = {p:p.read_bytes() for p in [rp, rb, mp, ap, lp]}; rb.unlink()
    audit_cases = ['raw_gram', 'projected_gram', 'raw_envelope', 'projected_envelope', 'diagnostic_rank',
        'retained_names', 'uniform_class', 'control_label', 'diagonal_hash', 'source_design_status',
        'promote_science', 'omit_last_audit', 'duplicate_audit', 'foreign_record_count']
    failure_evidence = 0
    for name in audit_cases:
        records = [json.loads(l) for l in gzip.decompress(originals[ap]).decode().splitlines()]
        r = next(r for r in records if r['numerical_audit'] is not None)
        value = r['numerical_audit']
        if name in ['raw_gram', 'projected_gram']: value[name][0][0] += .01
        elif name == 'raw_envelope': value['raw_roundoff_envelope'][0][0] *= 1.02
        elif name == 'projected_envelope': value['projected_roundoff_envelope'][0][0] *= 1.02
        elif name == 'diagnostic_rank': value['raw_diagnostics']['rank'] += 1
        elif name == 'retained_names': r['retained_kernel_names'] = list(reversed(r['retained_kernel_names']))
        elif name == 'uniform_class': r['residual_diagonal_is_exact_uniform_one'] = not r['residual_diagonal_is_exact_uniform_one']
        elif name == 'control_label': r['control_policy'] = 'foreign'
        elif name == 'diagonal_hash': r['diagonal_sha256'] = 'foreign'
        elif name == 'source_design_status': r['source_design_disposition'] = 'empty_setting'
        elif name == 'promote_science': r['scientific_eligibility'] = True
        elif name == 'omit_last_audit': records.pop()
        elif name == 'duplicate_audit': records[-1] = records[0]
        else: r['records'] += 1
        ap.write_bytes(gzip.compress((''.join(json.dumps(r) + '\n' for r in records)).encode(), mtime=0))
        m = json.loads(originals[mp]); receipt = json.loads(originals[rp])
        e = next(e for e in m if e['cohort_id'] == entry['cohort_id']); e['audit_sha256'] = sha(ap)
        write(mp, m); receipt['artifacts'][ap.relative_to(root).as_posix()] = sha(ap); receipt['artifacts'][mp.name] = sha(mp); write(rp, receipt)
        rejected(lambda: run(pp, True))
        failures = root / 'failures'
        if failures.exists():
            capture = a.output / ('preserved-rejected-' + name); failures.rename(capture)
            for f in capture.glob('*/failure.json'):
                data = json.loads(f.read_text()); assert data['scientific_eligibility'] is False
                assert (f.parent / 'original_numeric_inputs.npz').exists(); failure_evidence += 1
        for p in [rp, mp, ap]: p.write_bytes(originals[p])
    link_cases = ['foreign_link', 'missing_link', 'duplicate_link', 'changed_original_setting', 'changed_setting_ordinal']
    for name in link_cases:
        lines = gzip.decompress(originals[lp]).decode().splitlines(); fields = lines[0].split('\t')
        if name == 'missing_link': lines.pop()
        elif name == 'duplicate_link': lines[-1] = lines[1]
        else:
            cells = lines[1].split('\t')
            field = 'audit_id' if name == 'foreign_link' else 'scenario_id' if name == 'changed_original_setting' else 'source_setting_ordinal'
            cells[fields.index(field)] = 'foreign' if field != 'source_setting_ordinal' else str(int(cells[fields.index(field)]) + 1)
            lines[1] = '\t'.join(cells)
        lp.write_bytes(gzip.compress(('\n'.join(lines) + '\n').encode(), mtime=0))
        m = json.loads(originals[mp]); receipt = json.loads(originals[rp]); e = next(e for e in m if e['cohort_id'] == entry['cohort_id'])
        e['link_sha256'] = sha(lp); write(mp, m)
        receipt['artifacts'][lp.relative_to(root).as_posix()] = sha(lp); receipt['artifacts'][mp.name] = sha(mp); write(rp, receipt)
        rejected(lambda: run(pp, True))
        for p in [rp, mp, lp]: p.write_bytes(originals[p])
    replay = run(pp, True); previous = dict(reader); current = dict(replay)
    previous.pop('checked_utc'); current.pop('checked_utc'); assert previous == current
    verify(replay['source_hashes'])
    numerical = producer['numerical_audit_rows'] - producer['audit_status_counts'].get('insufficient_residual_dimension', 0)
    result = dict(status='passed_complete_declared_four_control_covariance_numerical_contracts_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), synthetic_cases=24, synthetic_cohorts=5,
        synthetic_designs=150, synthetic_settings=600, audit_records=6000, setting_audit_links=24000,
        numerical_calculations_completed=numerical, retained_basis_dimensions=[4, 5, 6],
        audit_status_counts=producer['audit_status_counts'], link_status_counts=producer['link_status_counts'],
        rehashed_audit_corruptions_rejected=audit_cases, rehashed_link_corruptions_rejected=link_cases,
        original_numeric_failure_captures=failure_evidence, completed_stage_restart_refusals=True,
        substantive_positive_restoration=True, source_and_journal_fixtures_synthetic=True,
        unchanged_comparison_rtol=3e-9, unchanged_comparison_atol=2e-8,
        source_hashes={**replay['source_hashes'], str(root / 'readback.json'):sha(root / 'readback.json')},
        scientific_eligibility=False, full_real_data_weighted_qualification_complete=False,
        scope=plan['scope'])
    with a.receipt.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_hashes'}), flush=True)


if __name__ == '__main__': main()
