#!/usr/bin/env python3
"""Complete declared source/export/readback contracts for positive diagonals."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from ancestral_chain_attempt import sha
from full_nonuniform_covariance_cone_sources import load
from full_retained_shared_entity_fit_fixtures_v2 import setup
from reference_measurement_union_sources import verify
from run_full_nonuniform_covariance_cones import run


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def rejected(action):
    try:
        action()
    except (AssertionError, ValueError, KeyError, StopIteration, FileExistsError):
        return
    raise AssertionError('Malformed full diagonal-cone proof accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args(); root = args.output.resolve(); root.mkdir(exist_ok=False)
    assert not args.receipt.exists()
    (root / 'fixture').mkdir()
    fp = setup(root / 'fixture'); fit = json.loads(fp.read_text())
    reduced = json.loads(Path(fit['retained_plan']).read_text())
    pp = root / 'plan.json'; output = root / 'cones'
    modules = ['nonuniform_covariance_cone', 'independent_nonuniform_covariance_cone',
               'full_nonuniform_covariance_cone_sources', 'run_full_nonuniform_covariance_cones',
               'check_full_nonuniform_covariance_cones_v2', 'full_retained_shared_entity_fit_fixtures_v2']
    plan = dict(parent_plan=reduced['exact_plan'], parent_completion=reduced['exact_completion'],
        expected=dict(logical_cases=24, cohorts=5, certificates=10), output=str(output),
        resources=dict(minimum_free_disk_gib=0),
        pins={'scripts/' + name + '.py': sha('scripts/' + name + '.py') for name in modules},
        scope='Complete declared five-cohort ten-certificate synthetic proof export/readback; '
              'parent source/journal fixtures synthetic, no real weighting or biological acceptance.')
    write(pp, plan)
    produced = run(pp); checked = run(pp, True)
    assert (produced['cohorts'], produced['certificates']) == (5, 10)
    assert produced['residual_diagonal_prepared'] is produced['raw_reml_basis_qualification_complete'] is False
    rejected(lambda: run(pp)); rejected(lambda: run(pp, True))
    data = output / 'cohort_cones.jsonl'; receipt = output / 'receipt.json'
    initial = {p: p.read_bytes() for p in [data, receipt, output / 'readback.json']}
    (output / 'readback.json').unlink()
    output_cases = ['omit_row', 'duplicate_row', 'parent_sha', 'forward', 'inverse', 'target_fold',
        'lost_kernel', 'negative_map', 'foreign_contract', 'promote_weighting', 'claim_diagonal', 'promote_science']
    for name in output_cases:
        rows = [json.loads(line) for line in initial[data].decode().splitlines()]
        m = rows[0]['variance_map']
        if name == 'omit_row': rows.pop()
        elif name == 'duplicate_row': rows[-1] = rows[0]
        elif name == 'parent_sha': rows[0]['parent_exact_certificate_sha256'] = 'foreign'
        elif name == 'forward': m['forward'][1][4] = .25
        elif name == 'inverse': m['nonnegative_right_inverse'][1][1] = 0.
        elif name == 'target_fold': m['forward'][0][1] = 1.
        elif name == 'lost_kernel': m['retained_names'].pop(1)
        elif name == 'negative_map': m['forward'][1][4] = -.5
        elif name == 'foreign_contract': rows[0]['source_contract'] = 'foreign'
        elif name == 'promote_weighting': m['nonuniform_weighting_accepted'] = True
        elif name == 'claim_diagonal': m['residual_diagonal_prepared'] = True
        else: rows[0]['scientific_eligibility'] = True
        data.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        value = json.loads(initial[receipt]); value['artifacts']['cohort_cones.jsonl'] = sha(data); write(receipt, value)
        rejected(lambda: run(pp, True))
        data.write_bytes(initial[data]); receipt.write_bytes(initial[receipt])
    cp = Path(plan['parent_completion']); completed = json.loads(cp.read_text())
    archive = Path(completed['full_hash_archive'])
    cert = Path(json.loads(Path(plan['parent_plan']).read_text())['output']) / 'cohort_certificates.jsonl'
    parents = {p: p.read_bytes() for p in [cp, archive, cert]}
    source_cases = ['foreign_completion', 'omit_certificate', 'duplicate_certificate',
                    'promote_parent', 'wrong_parent_mode', 'different_membership']
    for name in source_cases:
        closure = json.loads(parents[cp]); proof = json.loads(parents[archive])
        rows = [json.loads(line) for line in parents[cert].decode().splitlines()]
        if name == 'foreign_completion': closure['status'] = 'foreign'
        elif name == 'omit_certificate': rows.pop()
        elif name == 'duplicate_certificate': rows[-1] = rows[0]
        elif name == 'promote_parent': rows[0]['scientific_eligibility'] = True
        elif name == 'wrong_parent_mode': rows[0]['certificate']['loading_mode'] = 'foreign'
        else: rows[0]['cohort_rows_sha256'] = 'foreign'
        cert.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        proof['source_hashes'][str(cert)] = sha(cert); write(archive, proof)
        closure['full_hash_archive_sha256'] = sha(archive); write(cp, closure)
        rejected(lambda: load(plan, pp))
        for p, raw in parents.items():
            p.write_bytes(raw)
    replayed = run(pp, True)
    assert (output / 'readback.json').read_bytes() == initial[output / 'readback.json']
    verify(produced['source_hashes']); verify(replayed['source_hashes'])
    for p, raw in parents.items():
        assert p.read_bytes() == raw
    bindings = dict(replayed['source_hashes'])
    bindings.update({str(p): sha(p) for p in [data, receipt, output / 'readback.json', pp, Path(__file__)]})
    result = dict(status='passed_complete_declared_nonuniform_covariance_cone_source_export_contracts_v2',
        checked_utc=datetime.now(timezone.utc).isoformat(), synthetic_cohorts=5, synthetic_certificates=10,
        rehashed_output_cases_rejected=output_cases, rehashed_source_cases_rejected=source_cases,
        completed_producer_and_reader_restarts_rejected=True,
        byte_exact_source_and_positive_artifact_restoration=True, positive_full_bindings_rehashed=True,
        source_and_journal_fixtures_synthetic=True, source_hashes=bindings, scientific_eligibility=False,
        scope=plan['scope'])
    with args.receipt.open('x') as handle:
        json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'source_hashes'}))


if __name__ == '__main__':
    main()
