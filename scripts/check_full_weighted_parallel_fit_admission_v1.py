#!/usr/bin/env python3
"""Check complete parallel timing admission and private rebound corruptions."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil

from ancestral_chain_attempt import sha
from check_full_reduced_covariance_qualification import closed
from check_full_weighted_shared_entity_fits_v2 import rejected
from full_weighted_parallel_fit_admission_v1 import closed_timing, capacity, production_contract, FIELDS, COMPLETED
from prepare_weighted_timing_parallel_source_reference_v1 import READER
from prepare_baliphy_scalar_v6_preflight_v1 import project_sources
from reference_measurement_union_sources import bind, verify


def write(path, value):
    with Path(path).open('x') as handle:json.dump(value, handle, indent=2);handle.write('\n')


def overwrite(path, value):Path(path).write_text(json.dumps(value, indent=2)+'\n')


def copied(parent, destination):
    destination.mkdir(exist_ok=False)
    plan = json.loads(parent.read_text())
    original = Path(plan['output'])
    root = destination/'timing'
    shutil.copytree(original, root)
    tp = destination/'timing.plan.json'
    plan.update(output=str(root), scope='Private copied parallel timing custody fixture; '
        'actual native probe outputs inherited from the qualified parent, '
        'no producer/reader execution or original journal claim for rebound headers.')
    write(tp, plan)
    stage = json.loads((root/'stage_plan.json').read_text())
    stage['plan_sha256'] = sha(tp)
    overwrite(root/'stage_plan.json', stage)
    manifest = json.loads((root/'cohort_manifest.json').read_text())
    for index, part in enumerate(manifest):
        cp = root/part['receipt_path']
        value = json.loads(cp.read_text())
        value['stage'] = stage
        overwrite(cp, value)
        part['receipt_sha256'] = sha(cp)
        for role in ['producer', 'reader']:
            wp = root/'checkpoints'/(str(index).zfill(5)+'.'+role+'.json')
            worker = json.loads(wp.read_text())
            worker['manifest'] = part
            overwrite(wp, worker)
    overwrite(root/'cohort_manifest.json', manifest)
    producer = json.loads((root/'receipt.json').read_text())
    producer['plan_sha256'] = sha(tp)
    producer['source_hashes'].pop(str(parent))
    producer['source_hashes'][str(tp)] = sha(tp)
    producer['artifacts'] = {name:sha(root/name) for name in producer['artifacts']}
    overwrite(root/'receipt.json', producer)
    reader = json.loads((root/'readback.json').read_text())
    reader['plan_sha256'] = sha(tp)
    reader['producer_receipt_sha256'] = sha(root/'receipt.json')
    reader['source_hashes'] = dict(producer['source_hashes'])
    reader['source_hashes'][str(root/'receipt.json')] = sha(root/'receipt.json')
    reader['source_hashes'].update({str(root/name):digest for name,digest in producer['artifacts'].items()})
    reader['reader_checkpoint_artifacts'] = {name:sha(root/name) for name in reader['reader_checkpoint_artifacts']}
    overwrite(root/'readback.json', reader)
    overwrite(root/'reader_completed.json', dict(output=str(root/'readback.json'), sha256=sha(root/'readback.json'),
        status=READER, plan_sha256=sha(tp)))
    return tp


def completion(destination, tp, case=None):
    plan = json.loads(tp.read_text())
    root = Path(plan['output'])
    rp, rb = root/'receipt.json', root/'readback.json'
    producer, reader = [json.loads(path.read_text()) for path in [rp, rb]]
    fields = {key:reader[key] for key in FIELDS}
    paths = [tp, rp, rb, *[root/name for name in producer['artifacts']],
             *[root/name for name in reader['reader_checkpoint_artifacts']]]
    cp = closed(destination, 'closure', COMPLETED, paths, fields)
    value = json.loads(cp.read_text())
    value.update(producer_receipt=str(rp), producer_receipt_sha256=sha(rp),
                 independent_readback=str(rb), independent_readback_sha256=sha(rb))
    if case == 'unfinished_closure':value['exact_process_journals_checked'] = 1
    if case == 'wrong_completion_status':value['status'] = 'foreign'
    overwrite(cp, value)
    return cp


def corrupt(tp, case):
    plan = json.loads(tp.read_text())
    root = Path(plan['output'])
    rp, rb = root/'receipt.json', root/'readback.json'
    producer, reader = [json.loads(path.read_text()) for path in [rp, rb]]
    manifest = json.loads((root/'cohort_manifest.json').read_text())
    part = manifest[0]
    cp = root/part['receipt_path']
    fp, pp = root/part['census_path'], root/part['probes_path']
    if case == 'foreign_fit_contract':producer['fit_contract'] = reader['fit_contract'] = 'foreign'
    elif case in ['omit_candidate','duplicate_candidate','wrong_candidate_digest']:
        rows = [json.loads(line) for line in gzip.decompress(fp.read_bytes()).decode().splitlines()]
        if case == 'omit_candidate':rows.pop()
        elif case == 'duplicate_candidate':rows[-1] = deepcopy(rows[0])
        else:rows[0]['identity_sha256'] = 'foreign'
        fp.write_bytes(gzip.compress((''.join(json.dumps(row)+'\n' for row in rows)).encode(), mtime=0))
    elif case in ['foreign_representative','missing_group']:
        entry = next((p for p in manifest if json.loads((root/p['probes_path']).read_text())), None)
        if entry:
            part = entry
            cp = root/part['receipt_path']
            pp = root/part['probes_path']
            values = json.loads(pp.read_text())
            if case == 'missing_group':values.pop()
            else:values[0]['representative']['candidate_id'] = 'foreign'
            overwrite(pp, values)
        else:producer['timing_groups'] = reader['timing_groups'] = 1
    elif case == 'missing_numeric_replay':reader['fresh_numeric_probes_replayed'] += 1
    elif case == 'promoted_science':producer['scientific_eligibility'] = reader['scientific_eligibility'] = True
    elif case == 'altered_hardware_claim':reader['independent_hardware_timing_reimplementation'] = True
    elif case == 'changed_budget':producer['conditional_budget_weighted_seconds'] = reader['conditional_budget_weighted_seconds'] = 1.
    elif case == 'extra_artifact':
        extra = root/'invented.json'
        write(extra, dict(invented=True))
        producer['artifacts']['invented.json'] = sha(extra)
    elif case == 'wrong_stage_hash':producer['plan_sha256'] = reader['plan_sha256'] = 'foreign'
    elif case in ['unfinished_closure','wrong_completion_status','changed_reader_source']:pass
    elif case.startswith('checkpoint_'):
        wp = root/'checkpoints'/('00000.'+('producer' if case.startswith('checkpoint_producer') else 'reader')+'.json')
        value = json.loads(wp.read_text())
        if case.endswith('cache_false'):value['cached_numeric_inputs_preserved'] = False
        elif case.endswith('submitted_false'):value['submitted_numeric_inputs_preserved'] = False
        elif case.endswith('as_changed'):value['worker_address_space_limit_bytes'] += 1
        elif case.endswith('count_changed'):value['candidate_rows'] += 1
        elif case.endswith('index_changed'):value['cohort_index'] += 1
        elif case.endswith('pid_changed'):value['worker']['pid'] = 0
        elif case.endswith('cache_unchecked'):value['cached_input_array_bindings_checked'] = 0
        elif case.endswith('task_unchecked'):value['submitted_task_array_bindings_checked'] = 0
        elif case.endswith('science_true'):value['scientific_eligibility'] = True
        elif case.endswith('review_invented'):value['review_artifacts']['invented'] = 'foreign'
        else:raise AssertionError(case)
        overwrite(wp, value)
    elif case == 'missing_reader_checkpoint_binding':reader['reader_checkpoint_artifacts'].pop('checkpoints/00000.reader.json')
    else:raise AssertionError(case)
    # Rehash custody after mutations; semantic guards must still reject.
    value = json.loads(cp.read_text())
    value.update(census_sha256=sha(root/part['census_path']), probes_sha256=sha(root/part['probes_path']))
    overwrite(cp, value)
    for key in ['census','probes','receipt']:part[key+'_sha256'] = sha(root/part[key+'_path'])
    overwrite(root/'cohort_manifest.json', manifest)
    index = next(i for i,p in enumerate(manifest) if p['cohort_id'] == part['cohort_id'])
    for role in ['producer','reader']:
        wp = root/'checkpoints'/(str(index).zfill(5)+'.'+role+'.json')
        value = json.loads(wp.read_text())
        value['manifest'] = part
        overwrite(wp, value)
    producer['artifacts'] = {name:sha(root/name) for name in producer['artifacts']}
    overwrite(rp, producer)
    reader['producer_receipt_sha256'] = sha(rp)
    reader['source_hashes'] = dict(producer['source_hashes'])
    reader['source_hashes'][str(rp)] = sha(rp)
    reader['source_hashes'].update({str(root/name):digest for name,digest in producer['artifacts'].items()})
    reader['reader_checkpoint_artifacts'] = {name:sha(root/name) for name in reader['reader_checkpoint_artifacts']}
    if case == 'changed_reader_source':reader['source_hashes'][next(iter(reader['source_hashes']))] = 'foreign'
    overwrite(rb, reader)
    overwrite(root/'reader_completed.json', dict(output=str(rb), sha256=sha(rb), status=READER, plan_sha256=sha(tp)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(exist_ok=False)
    assert not a.receipt.exists()
    gate = Path('metadata/full_weighted_timing_parallel_software_validation_20261004_v1.json')
    qualified = json.loads(gate.read_text())
    verify(qualified['source_hashes'])
    pins = {}
    project_sources(pins, [Path(__file__)])
    bind(pins, gate)
    fixture = Path('data/software_audits/full-weighted-timing-parallel-20261004-v1')
    names = ['unqualified','qualified-common-pair','qualified-pair-exception']
    bad = []
    grids = []
    cases = ['unfinished_closure','wrong_completion_status','foreign_fit_contract','changed_reader_source',
        'omit_candidate','duplicate_candidate','wrong_candidate_digest','foreign_representative','missing_group',
        'missing_numeric_replay','promoted_science','altered_hardware_claim','changed_budget','extra_artifact','wrong_stage_hash',
        'checkpoint_producer_cache_false','checkpoint_producer_submitted_false','checkpoint_reader_cache_false',
        'checkpoint_reader_submitted_false','checkpoint_reader_as_changed','checkpoint_reader_count_changed',
        'checkpoint_reader_index_changed','checkpoint_reader_pid_changed','checkpoint_reader_cache_unchecked',
        'checkpoint_reader_task_unchecked','checkpoint_reader_science_true','checkpoint_reader_review_invented',
        'missing_reader_checkpoint_binding']
    for name in names:
        parent = fixture/(name+'-parallel.timing.plan.json')
        folder = a.output/(name+'-positive')
        tp = copied(parent, folder)
        cp = completion(folder, tp)
        plan = json.loads(tp.read_text())
        source, hashes, result = closed_timing(plan['fit_plan'], tp, cp)
        assert result['fresh_admission_candidate_rows'] == 24000 and result['full_parallel_checkpoint_pairs_checked'] == 5
        assert not result['admission_repeats_numeric_probes'] and not result['admission_retains_all_probe_arrays_in_memory']
        fit = json.loads(Path(plan['fit_plan']).read_text())
        fit['resources'] = dict(minimum_free_disk_gib=3172, output_scratch_reserve_gib=3072)
        resources = capacity(result, fit)
        assert resources['resources_installed'] is False and resources['production_finish_eta'] is None
        rejected(lambda:production_contract(fit, result))
        grids.append(dict(name=name, candidate_rows=24000, timing_groups=result['timing_groups'], checked_checkpoint_pairs=5))
        for path, digest in hashes.items():bind(pins, path, digest)
        for case in cases:
            destination = a.output/(name+'-private-'+case)
            tp = copied(parent, destination)
            corrupt(tp, case)
            cp = completion(destination, tp, case)
            plan = json.loads(tp.read_text())
            rejected(lambda:closed_timing(plan['fit_plan'], tp, cp))
            bad.append([name, case])
    full = Path('metadata/full_weighted_timing_parallel_plan_20261004_v1.json')
    real = json.loads(full.read_text())
    assert not Path(real['completion']).exists()
    rejected(lambda:closed_timing(real['fit_plan'], full, real['completion']))
    assert not Path(json.loads(Path(real['fit_plan']).read_text())['output']).exists()
    bind(pins, full)
    for path in a.output.rglob('*'):
        if path.is_file():bind(pins, path)
    verify(qualified['source_hashes'])
    verify(pins)
    proof = dict(status='passed_complete_parallel_timing_checkpoint_fit_admission_v1',
        checked_utc=datetime.now(timezone.utc).isoformat(), complete_synthetic_grids=grids,
        total_admission_candidate_rows=72000, full_checkpoint_pairs_checked=15,
        private_rehashed_corruptions_rejected=bad, corruption_cases_per_grid=len(cases),
        unchanged_qualified_parent_bytes_verified=True, positive_custody_copies_passed=True,
        actual_missing_timing_prerequisite_refused=True, synthetic_scope_cannot_admit_production=True,
        full_native_timing_probes_repeated=False, full_probe_arrays_accumulated=False,
        source_and_journal_closure_fixtures_synthetic=True, source_hashes=pins,
        production_fitting_launched=False, scientific_eligibility=False, gpu=False,
        scope='All three complete parallel timing grids/72000candidate identities and15producer/reader '
              'checkpoint pairs reconstructed. Actual parent320native groups already qualified; '
              'no fresh timing/optimizer run. Eighty-four distinct private copied rehashed '
              'accounting/checkpoint/closure corruptions reject without changing parent bytes. '
              'Cohort-wise planning reduction preserves counts/costs without retaining all '
              'probe arrays. Synthetic closures do not satisfy production scope; real missing '
              'timing refuses fitting. Capacity allocations remain neither ETA nor inference.')
    write(a.receipt, proof)
    print(json.dumps({k:v for k,v in proof.items() if k not in ['source_hashes','scope']}, indent=2))


if __name__ == '__main__':main()
