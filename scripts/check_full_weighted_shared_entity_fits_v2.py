#!/usr/bin/env python3
"""Complete weighted export/checkpoint/link contracts and actual-fit replays.

Native fits in the complete export grids are explicitly mocked; a separate
64-case serialized native-fit replay uses the frozen actual-D numeric reader.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from ancestral_chain_attempt import sha
from check_full_shared_entity_fits_v2 import settings
from check_retained_shared_entity_candidates import inputs
from check_weighted_shared_entity_candidates import source_grids
from full_weighted_fit_exports import atomic, validate_export
from prepare_full_weighted_shared_entity_fits import run
from readback_full_weighted_shared_entity_fits import run as readback
from readback_weighted_shared_entity_candidate import numeric
from reference_measurement_union_sources import verify
from weighted_shared_entity_candidate import READY, validate_source


def rejected(action):
    try: action()
    except (AssertionError, ValueError, ArithmeticError, KeyError, FileExistsError, FileNotFoundError, StopIteration): return
    raise AssertionError('Malformed weighted fitting evidence accepted')


def synthetic_candidate(source, plan, rows, identity, x, y, operators, audit, route, diagonal):
    validate_source(identity, audit, route, diagonal)
    if identity['source_combined_disposition'] != READY:
        return dict(**identity, disposition=identity['source_combined_disposition'], fit=None, numerical_attempted=False)
    return dict(**identity, disposition='shared_entity_fit_error_requires_review', fit=None, numerical_attempted=True,
        error_type='ArithmeticError', error_message='Explicit synthetic mocked native failure; not a scientific fit')


def synthetic_reader(source, plan, rows, identity, x, y, exported, operators, audit, route, diagonal):
    validate_source(identity, audit, route, diagonal); validate_export(identity, exported)
    if identity['source_combined_disposition'] != READY:
        return numeric(source, plan, rows, identity, x, y, exported, operators, audit, route, diagonal)
    assert exported == synthetic_candidate(source, plan, rows, identity, x, y, operators, audit, route, diagonal)
    return dict(disposition='original_numerical_failure_reproduced_requires_review', scientific_eligibility=False)


def plan_for(pp, output):
    value = json.loads(pp.read_text()); value.update(output=str(output), resources=dict(minimum_free_disk_gib=0),
        scope='Complete synthetic source/candidate/link/checkpoint grid; parent closures/journals synthetic, native producer/reader explicitly mocked. No production biological fit or pilot.')
    pp.write_text(json.dumps(value, indent=2) + '\n')


def corruptions(pp, root):
    receipt_path = root / 'receipt.json'; manifest_path = root / 'cohort_manifest.json'
    receipt_bytes = receipt_path.read_bytes(); manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes); first = manifest[0]
    fp = root / first['candidates']; lp = root / first['links']; cp = root / first['receipt']
    originals = {p:p.read_bytes() for p in [fp, lp, cp, receipt_path, manifest_path]}
    rejections = []; captures = 0
    for kind in ['omit_candidate','duplicate_candidate','foreign_candidate','changed_diagonal','promote_science','false_fit_failure',
                 'missing_link','duplicate_link','foreign_link','changed_ordinal','changed_original_setting','changed_control',
                 'missing_artifact','changed_source_hash','extra_artifact','false_stage_contract']:
        candidate_rows = [json.loads(l) for l in gzip.decompress(originals[fp]).decode().splitlines()]
        link_lines = gzip.decompress(originals[lp]).decode().splitlines()
        receipt = json.loads(receipt_bytes); current_manifest = json.loads(manifest_bytes); checkpoint = json.loads(originals[cp])
        if kind == 'omit_candidate': candidate_rows.pop()
        elif kind == 'duplicate_candidate': candidate_rows[-1] = deepcopy(candidate_rows[0])
        elif kind == 'foreign_candidate': candidate_rows[0]['candidate_id'] = 'foreign'
        elif kind == 'changed_diagonal': candidate_rows[0]['diagonal_sha256'] = 'foreign'
        elif kind == 'promote_science': candidate_rows[0]['scientific_eligibility'] = True
        elif kind == 'false_fit_failure':
            candidate_rows[0].update(disposition='shared_entity_fit_error_requires_review', numerical_attempted=True,
                error_type='ArithmeticError', error_message='invented')
        elif kind == 'missing_link': link_lines.pop()
        elif kind == 'duplicate_link': link_lines[-1] = link_lines[1]
        elif kind in ['foreign_link','changed_ordinal','changed_original_setting','changed_control']:
            fields = link_lines[0].split('\t'); cells = link_lines[1].split('\t')
            field = {'foreign_link':'candidate_id','changed_ordinal':'source_setting_ordinal',
                'changed_original_setting':'scenario_id','changed_control':'control_policy'}[kind]
            cells[fields.index(field)] = 'foreign'; link_lines[1] = '\t'.join(cells)
        elif kind == 'missing_artifact': receipt['artifacts'].pop(first['links'])
        elif kind == 'changed_source_hash': receipt['source_hashes'][next(iter(receipt['source_hashes']))] = 'foreign'
        elif kind == 'extra_artifact': receipt['artifacts']['invented'] = 'foreign'
        else: receipt['source_contract'] = 'foreign'
        fp.write_bytes(gzip.compress((''.join(json.dumps(r) + '\n' for r in candidate_rows)).encode(), mtime=0))
        lp.write_bytes(gzip.compress(('\n'.join(link_lines) + '\n').encode(), mtime=0))
        checkpoint.update(candidate_sha256=sha(fp), links_sha256=sha(lp))
        cp.write_text(json.dumps(checkpoint) + '\n')
        current_manifest[0].update(candidates_sha256=sha(fp), links_sha256=sha(lp), receipt_sha256=sha(cp))
        manifest_path.write_text(json.dumps(current_manifest) + '\n')
        for p in [fp, lp, cp, manifest_path]:
            if p.relative_to(root).as_posix() in receipt['artifacts']: receipt['artifacts'][p.relative_to(root).as_posix()] = sha(p)
        receipt_path.write_text(json.dumps(receipt) + '\n')
        rejected(lambda: readback(pp, root / 'rejected-reader.json'))
        independent = root / 'independent'
        if independent.exists():
            failures = list((independent / 'failures').glob('*/failure.json'))
            for f in failures:
                assert json.loads(f.read_text())['scientific_eligibility'] is False
                with np.load(f.parent / 'original_numeric_inputs.npz', allow_pickle=False) as values:
                    assert set(values.files) == {'case_rows','diagonal','design','response','species_factor'}
                captures += 1
            independent.rename(root / ('preserved-rejected-' + kind))
        for p, b in originals.items(): p.write_bytes(b)
        rejections.append(kind)
    return rejections, captures


def actual_serialized(root):
    gate_path = Path('metadata/weighted_shared_entity_candidate_software_validation_20261004_v1.json')
    gate = json.loads(gate_path.read_text()); verify(gate['source_hashes'])
    path = Path('data/software_audits/weighted-shared-entity-candidates-20261004-v1/actual_weighted_candidate_cases.json')
    values = json.loads(path.read_text()); assert len(values) == 64
    plan = settings(); results = []
    for row in values:
        source, rows, x, y, ops, retained, old, cert = inputs(row['pair_exception'], row['loading_mode'])
        from scipy import sparse
        operators = {} if row['policy'] == 'uniform' else {'target_node':sparse.eye(len(rows), format='csr')}
        operators.update(ops)
        response = y if row['outcome'] == 'rmsd_delta' else np.tanh(y) + .15*x[:,1]
        diagonal = np.asarray(row['diagonal']); value = json.loads(json.dumps(row['candidate'], allow_nan=False))
        identity = {k:v for k,v in value.items() if k not in ['disposition','fit','numerical_attempted','error_type','error_message']}
        audit = row['audit']; route = dict(names=audit['retained_kernel_names'],
            certificate_sha256=audit['exact_certificate_sha256'], route=audit['basis_route'],
            exact_uniform_one=audit['residual_diagonal_is_exact_uniform_one'])
        validate_export(identity, value)
        checked = numeric(source, plan, rows, identity, x, response, value, operators, audit, route, diagonal)
        assert checked == row['independent']
        results.append(dict(candidate_id=identity['candidate_id'], independent=checked))
    artifact = root / 'actual_serialized_fit_replays.json'; atomic(artifact, results)
    return {str(gate_path):sha(gate_path), str(path):sha(path), str(artifact):sha(artifact)}, results


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True); a = p.parse_args(); a.output.mkdir(exist_ok=False); assert not a.receipt.exists()
    setup = a.output / 'source-fixtures'; setup.mkdir()
    grids, bindings, source_rejections = source_grids(setup)
    completed = []; rejections = []; captures = 0
    with patch('prepare_full_weighted_shared_entity_fits.candidate', side_effect=synthetic_candidate), \
         patch('readback_full_weighted_shared_entity_fits.numeric', side_effect=synthetic_reader):
        for number in range(3):
            pp = setup / ('fit-plan-' + str(number) + '.json'); root = a.output / ('exports-' + str(number)); plan_for(pp, root)
            if number == 0:
                try: run(pp, stop_after_cohorts=2)
                except InterruptedError: pass
                else: raise AssertionError('Producer interruption not exercised')
                previous = {str(f):sha(f) for f in (root / 'cohorts').glob('*') if f.is_file()}
                assert len(list((root / 'cohorts').glob('*.receipt.json'))) == 2
            producer = run(pp)
            if number == 0: assert all(sha(f) == h for f,h in previous.items())
            assert producer['candidate_rows'] == 24000 and producer['setting_fit_links'] == 48000
            if number == 0:
                bad, n = corruptions(pp, root); rejections.extend(bad); captures += n
                try: readback(pp, root / 'readback.json', stop_after_cohorts=2)
                except InterruptedError: pass
                else: raise AssertionError('Reader interruption not exercised')
                prior = {str(f):sha(f) for f in (root / 'independent').glob('*') if f.is_file()}
                assert len(list((root / 'independent').glob('*.receipt.json'))) == 2
                # A newly hashed reader checkpoint is still replayed, and
                # cannot promote an original review or invent an accepted fit.
                ap = next((root / 'independent').glob('*.audits.jsonl.gz'))
                ip = ap.with_name(ap.name.replace('.audits.jsonl.gz', '.receipt.json'))
                before_ap, before_ip = ap.read_bytes(), ip.read_bytes()
                audit_rows = [json.loads(line) for line in gzip.decompress(before_ap).decode().splitlines()]
                audit_rows[0]['disposition'] = 'independently_audited_working_candidate_pending_inferential_calibration'
                ap.write_bytes(gzip.compress((''.join(json.dumps(r)+'\n' for r in audit_rows)).encode(),mtime=0))
                saved = json.loads(before_ip); saved['audits_sha256'] = sha(ap); ip.write_text(json.dumps(saved)+'\n')
                rejected(lambda: readback(pp, root / 'rejected-reader-checkpoint.json'))
                failures = root / 'independent' / 'failures'
                assert len(list(failures.glob('*/failure.json'))) == 1
                failures.rename(root / 'preserved-reader-checkpoint-rejection')
                ap.write_bytes(before_ap); ip.write_bytes(before_ip)
                rejections.append('rehashed_reader_checkpoint_promotes_review'); captures += 1
            reader = readback(pp, root / 'readback.json')
            if number == 0: assert all(sha(f) == h for f,h in prior.items())
            rejected(lambda: run(pp)); rejected(lambda: readback(pp, root / 'other-reader.json'))
            verify(reader['source_hashes']); bindings.update(reader['source_hashes'])
            bindings[str(root / 'readback.json')] = sha(root / 'readback.json')
            completed.append(dict(grid=number,candidate_rows=24000,setting_fit_links=48000,
                producer_status_counts=producer['candidate_status_counts'],independent_status_counts=reader['independent_candidate_status_counts'],
                producer=str(root / 'receipt.json'), reader=str(root / 'readback.json')))
    actual_bindings, replays = actual_serialized(a.output); bindings.update(actual_bindings)
    modules = ['check_full_weighted_shared_entity_fits_v2','check_full_weighted_shared_entity_fits','full_weighted_fit_exports','prepare_full_weighted_shared_entity_fits',
        'readback_full_weighted_shared_entity_fits','full_weighted_shared_entity_fit_sources','weighted_shared_entity_candidate',
        'readback_weighted_shared_entity_candidate','check_weighted_shared_entity_candidates','reference_measurement_union_sources']
    bindings.update({'scripts/'+m+'.py':sha('scripts/'+m+'.py') for m in modules}); verify(bindings)
    result = dict(status='passed_complete_four_control_fit_exports_readback_and_checkpoint_contracts_v2',
        checked_utc=datetime.now(timezone.utc).isoformat(),complete_synthetic_grids=completed,total_candidate_rows=72000,
        total_setting_fit_links=144000,actual_serialized_numeric_replays=64,rehashed_output_alterations_rejected=rejections,
        synthetic_source_closure_alterations_rejected=source_rejections,exact_numeric_failure_captures=captures,
        producer_closed_chunks_reused_without_rewrite=True,reader_closed_chunks_replayed_and_reused_without_rewrite=True,
        completed_producer_and_reader_restarts_refused=True,source_and_journal_fixtures_synthetic=True,
        complete_export_grids_native_producer_and_reader_explicitly_mocked=True,
        actual_serialized_numeric_replays_use_frozen_independent_reader=True,
        source_hashes=bindings,scientific_eligibility=False,production_fitting_launched=False,
        scope='Complete three declared five-cohort original four-control export/link/checkpoint grids,72000candidate rows/144000links, with explicitly mocked native failure/reader branches. Separate64actual saved numerical fits replay through frozen independent actual-D spectral/start/curvature/search reader. Every original setting ordinal retained. Closed chunks checked/replayed; failures and partial artifacts preserved. No real source qualification, timing, production fit, biological pilot or accepted effect.')
    atomic(a.receipt,result); print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','scope']},indent=2),flush=True)


if __name__ == '__main__': main()
