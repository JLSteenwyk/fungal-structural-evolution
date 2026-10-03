#!/usr/bin/env python3
"""Queue full native-output validation after the original sampler closes."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import sha
from independent_native_ancestral_alignment import fasta_records
from independent_short_sampler_outputs import NATIVE_ALPHABET, SUMMARY_FIELDS
from launch_baliphy_reference_sampler_qualification import launch
from launch_full_triad_sequence_geometry_followup import create
from record_project_runtime_checkpoint_v4 import fingerprint, journal_terminal
from run_baliphy_reference_preflight import verify
from run_independent_short_sampler_replay import PRODUCER_STATUS, READER_STATUS, COMPLETED_STATUS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', type=Path, required=True)
    args = parser.parse_args(); gate = json.loads(args.validation.read_text())
    assert gate['status'] == 'passed_full_independent_short_sampler_replay_software_contracts'
    assert gate['full_real_roles'] == 1620 and gate['full_real_groups'] == 405
    assert gate['full_producer_and_reader_serialization_checked']
    assert gate['actual_existing_native_roles'] == 3 and gate['actual_existing_native_alignments'] == 9
    assert gate['native_alphabet_verified'] == NATIVE_ALPHABET
    assert gate['artificial_failed_roles_retained'] == 2 and gate['intact_artificial_roles_in_unresolved_groups'] == 6
    assert len(gate['malformed_property_cases_rejected']) == 15 and len(gate['serialization_cases_rejected']) == 10
    verify({'pins': gate['source_hashes']})
    native_path = Path('metadata/baliphy_reference_sampler_qualification_plan_20261003.json')
    native = json.loads(native_path.read_text()); verify(native)
    dependency = native['launch_inventory']; inventory = json.loads(Path(dependency).read_text())
    original_closer = inventory['launches'][-1]; original = json.loads(Path(original_closer).read_text())
    original['launch'] = original_closer
    if fingerprint(original) is None: journal_terminal(original)
    assert psutil.virtual_memory().available >= 16 * 2**30 and shutil.disk_usage('.').free >= 64 * 2**30
    jobs = json.loads(Path(native['jobs']).read_text()); residues = {}
    for job in jobs:
        path = job['chain']['alignment']
        if path not in residues:
            with Path(path).open() as handle:
                residues[path] = sum(len(s.replace('-', '')) for s in fasta_records(handle).values())
    potential = sum(3 * 4 * residues[j['chain']['alignment']] for j in jobs)
    maximum = max(3 * 4 * residues[j['chain']['alignment']] for j in jobs)
    own = ['independent_short_sampler_outputs', 'run_independent_short_sampler_replay',
           'check_independent_short_sampler_replay', 'prepare_and_launch_independent_short_sampler_replay',
           'record_independent_short_sampler_replay_checkpoint', 'independent_native_ancestral_alignment',
           'independent_native_ancestral_topology']
    paths = [Path('scripts') / (name + '.py') for name in own]
    paths += [args.validation, native_path, Path(native['jobs']), Path(native['mapping']), Path(original_closer)]
    native_api = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/lib/bali-phy/haskell')
    paths += [native_api / name for name in ['Bio/Alphabet.hs', 'Bio/Alignment.hs', 'SModel/ASRV.hs',
              'Probability/Logger.hs', 'Probability/Distribution/PhyloCTMC/VariableA/Properties.hs']]
    pins = dict(native['pins']); pins.update({str(p): sha(p) for p in paths})
    output = 'results/ancestral/full-independent-short-sampler-replay-20261003-v1'
    scope = ('All1620originalshortsamplerroles,405quartets,135effectiveinputs,324aliases and every failed/invalid '
        'original disposition retained after complete original sampler source/artifact/two-journal closure. '
        'Independent manual source/runtime clades, candidate mappings, saved FASTA decoding and tip-residue anchored '
        'projections for all0/10/20frames; exact full serialized/NPZ replay. Every available native tip category/state '
        'array,4x20rate property matrix and empty condition schema checked. Internal category labels are not logged, '
        'are explicitly unavailable and never imputed from tips. No complete ancestral category trajectory, '
        'gamma-discretization likelihood proof, convergence, allocation repair, biological root/model/predictor '
        'acceptance, new inference/GPU/costs or scientific aim completion. Reader shares independent decoder; '
        'complete source/artifact/two original journals required for output accounting closure.')
    plan_path = 'metadata/independent_short_sampler_replay_plan_20261003.json'
    plan = dict(sampler_plan=str(native_path), mapping=native['mapping'], pins=pins, output=output,
        completion='metadata/independent_short_sampler_replay_completed_20261003.json',
        launch_inventory='metadata/independent_short_sampler_replay_launches_20261003.json',
        resources=dict(checked_utc=datetime.now(timezone.utc).isoformat(), cpus=2, memory_gib=16, swap_gib=0,
            blas_threads=1, maximum_saved_alignments=4860, maximum_candidate_frames=19440,
            potential_projected_uint8_bytes=potential, maximum_one_role_projected_uint8_bytes=maximum,
            planning_output_gib=8, minimum_free_disk_gib=64, active_hours_per_stage=[.05,8], runtime_uncalibrated=True,
            finish_eta=None, new_native_inference_jobs=0, gpu=False, new_cost_usd=0,
            available_memory_gib=psutil.virtual_memory().available/2**30, available_disk_gib=shutil.disk_usage('.').free/2**30,
            caveat='Exact potential projection bytes assume all roles pass; native JSON/FASTA workspace and alignment expansion remain unbounded by this estimate. Enforced cgroup16GiB/no swap and free-disk guard are distinct from output/runtime planning. Each stage rereads closed original artifacts; no ETA or posterior claim.'), scope=scope)
    assert not Path(output).exists(); create(plan_path, plan)
    producer = launch('independent-short-sampler-replay',
        [sys.executable, 'scripts/run_independent_short_sampler_replay.py', '--plan', plan_path],
        [original_closer], plan_path, cpus=2, memory=16)
    reader = launch('independent-short-sampler-replay-readback',
        [sys.executable, 'scripts/run_independent_short_sampler_replay.py', '--plan', plan_path, '--reader'],
        [producer], plan_path, cpus=2, memory=16)
    completion_path = 'metadata/independent_short_sampler_replay_completion_plan_20261003.json'
    spec = dict(source_plan=plan_path, producer_receipt=output+'/receipt.json', independent_readback=output+'/readback.json',
        producer_status=PRODUCER_STATUS, reader_status=READER_STATUS, completed_status=COMPLETED_STATUS,
        summary_fields=SUMMARY_FIELDS, launches=[producer,reader],
        pins={p:sha(p) for p in [plan_path,producer,reader,'scripts/close_full_triad_sequence_stage.py',
              'scripts/record_completed_process_handoffs_v2.py']}, output=plan['completion'], scope=scope)
    create(completion_path,spec)
    closer = launch('independent-short-sampler-replay-closure',
        [sys.executable,'scripts/close_full_triad_sequence_stage.py','--plan',completion_path],
        [producer,reader],completion_path,cpus=2,memory=16)
    create(plan['launch_inventory'],dict(status='queued_full_independent_short_sampler_native_replay',
        source_plan=plan_path,source_plan_sha256=sha(plan_path),launches=[producer,reader,closer],
        dependency_original_sampler_closure=original_closer,new_native_inference_jobs=0,gpu=False,new_cost_usd=0))


if __name__ == '__main__': main()
