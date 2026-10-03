#!/usr/bin/env python3
"""Prepare all 1,620 startup roles and immutable programs; no posterior launch."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

from Bio import Phylo
import psutil

from ancestral_chain_attempt import sha, write_json
from baliphy_reference_initialization import alignment_records, transform


SCOPE = ('Full135 effective inputs x3 original priors x4 fresh seed roles, all324 aliases. '
         'Native --test startup only: supplied alignment is a starting state, not observed aligned data. '
         'Original distribution, density, kernels, fixed-tip/free-ancestor representation and priors retained. '
         'All native failures and invalid startups retained; no posterior sampling, memory repair, mixing '
         'qualification, GPU, new costs or scientific acceptance. Complete serialized/source/journal closure required.')


def build(chains, proposals, out, binary, prlimit, api_paths):
    assert len(chains) == len({x['chain_id'] for x in chains}) == 1620
    proposed = {x['source_chain_id']: x for x in proposals}
    assert len(proposed) == len(proposals) == 1620 and set(proposed) == {x['chain_id'] for x in chains}
    assert len({x['fresh_seed'] for x in proposals}) == 1620
    assert all(isinstance(x['fresh_seed'], int) and 1 <= x['fresh_seed'] <= 2**31 - 1 for x in proposals)
    assert {x['fresh_seed'] for x in proposals}.isdisjoint({x['seed'] for x in chains})
    assert len({x['effective_input_group'] for x in chains}) == 135
    assert Counter(x['prior_label'] for x in chains) == {'broad': 540, 'centered': 540, 'package': 540}
    assert len({alias for x in chains for alias in x['original_configuration_ids']}) == 324
    groups = defaultdict(list); originals = {}; jobs = []; bindings = {}; matrices = {}; trees = {}
    for chain in chains: groups[chain['effective_input_group'] + '-' + chain['prior_label']].append(chain)
    assert len(groups) == 405
    for group, quartet in sorted(groups.items()):
        assert len(quartet) == 4 and {x['chain'] for x in quartet} == {1, 2, 3, 4}
        assert len({x['program'] for x in quartet}) == 1
        shared = ['effective_input_group', 'prior_label', 'original_configuration_ids',
                  'family', 'proteins', 'alignment', 'alignment_sha256', 'tree',
                  'tree_sha256', 'program_sha256']
        assert all(all(x[k] == quartet[0][k] for k in shared) for x in quartet)
        chain = quartet[0]; original = Path(chain['program'])
        assert sha(original) == chain['program_sha256']
        assert sha(chain['alignment']) == chain['alignment_sha256']
        assert sha(chain['tree']) == chain['tree_sha256']
        for name in ['program', 'alignment', 'tree']:
            bindings[str(Path(chain[name]))] = chain[name + '_sha256']
        if chain['alignment'] not in matrices: matrices[chain['alignment']] = alignment_records(chain['alignment'])
        reference = matrices[chain['alignment']]
        if chain['tree'] not in trees: trees[chain['tree']] = Phylo.read(chain['tree'], 'newick')
        tree = trees[chain['tree']]; tips = {x.name for x in tree.get_terminals()}
        assert tips == {x[0] for x in reference} and len(tips) == chain['proteins']
        program = out / 'models' / (group + '.hs'); program.parent.mkdir(exist_ok=True)
        with program.open('x') as f: f.write(transform(original.read_text()))
        originals[group] = dict(source_program=str(original), source_program_sha256=sha(original),
                                reference_program=str(program), reference_program_sha256=sha(program),
                                reversible_initialization_edits_only=True)
        for chain in sorted(quartet, key=lambda x: x['chain']):
            q = proposed[chain['chain_id']]
            assert q['source_seed'] == chain['seed'] and q['model_input_identity'] == group
            assert q['original_configuration_ids'] == chain['original_configuration_ids']
            assert q['iterations'] == 10000 and q['production_launch_allowed'] is False
            paths = [binary, prlimit, program, Path(chain['alignment']).resolve(), Path(chain['tree']).resolve(), *api_paths]
            config = dict(command=[str(prlimit), '--as=' + str(12 * 2**30), '--cpu=600',
                '--fsize=' + str(256 * 2**20), '--', str(binary), '--seed', str(q['fresh_seed']),
                'run', str(program), '--test', '--log-format', 'json'],
                timeout_seconds=900, pins={str(p.resolve()): sha(p) for p in paths})
            jobs.append(dict(chain=chain, fresh_seed=q['fresh_seed'], program=str(program), config=config))
    assert len(jobs) == 1620 and len(originals) == 405
    return jobs, originals, bindings


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--software-validation', type=Path, required=True)
    p.add_argument('--grid-validation', type=Path, required=True)
    a = p.parse_args()
    validation = json.loads(a.software_validation.read_text())
    assert validation['status'] == 'passed_native_reference_initialization_software_contracts'
    assert (validation['native_initializations'], validation['original_initializations'], validation['kernel_fixtures']) == (12, 12, 3)
    assert len(validation['invalid_inputs_and_audits_rejected']) == 23
    for path, h in validation['source_hashes'].items(): assert sha(path) == h, path
    grid = json.loads(a.grid_validation.read_text())
    assert grid['status'] == 'passed_full_reference_startup_grid_software_contracts'
    assert grid['full_chains'] == 1620 and grid['full_quartets'] == 405
    assert len(grid['altered_designs_and_missing_dispositions_rejected']) == 11
    for path, h in grid['source_hashes'].items(): assert sha(path) == h, path
    chain_path = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    closure_path = Path('metadata/baliphy_horizon_resources_completed_20261003.json')
    closed = json.loads(closure_path.read_text())
    assert closed['status'] == 'complete_verified_full_baliphy_horizon_resource_inventory'
    archive = Path(closed['full_hash_archive']); assert sha(archive) == closed['full_hash_archive_sha256']
    archived = json.loads(archive.read_text())
    proposals_path = archive.parent / 'proposed.jsonl'
    assert sha(proposals_path) == archived['source_hashes'][str(proposals_path)]
    source_manifest = Path('data/software_audits/baliphy-4.3-reference-initialization-20261003/source_manifest.json')
    sources = json.loads(source_manifest.read_text())
    api_paths = []
    for name, r in sources['files'].items():
        installed = Path(r['installed_path']); downloaded = source_manifest.parent / name
        assert sha(installed) == sha(downloaded) == r['sha256']
        api_paths.append(installed.resolve())
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    assert sha(binary) == validation['binary_sha256']
    jobs, models, bindings = build(json.loads(chain_path.read_text()),
        [json.loads(x) for x in proposals_path.read_text().splitlines()], out, binary, Path('/usr/bin/prlimit'), api_paths)
    write_json(out / 'jobs.json', jobs); write_json(out / 'model_transformations.json', models)
    old = Path('results/ancestral/baliphy-prior-initialization-20260927-v1')
    earlier = [json.loads(p.read_text()) for p in old.glob('*/receipt.json')]
    assert len(earlier) == 405 and all(x['exit_code'] == 0 for x in earlier)
    measured = sum(x['elapsed_seconds'] for x in earlier) * 4
    resources = dict(checked_utc=datetime.now(timezone.utc).isoformat(), workers=2, cpus=2,
        memory_gib=24, swap_gib=0, blas_threads=1, per_worker_address_space_gib=12,
        per_worker_cpu_seconds_cap=600, per_worker_wall_seconds_cap=900,
        per_file_limit_mib=256, output_allowance_gib=16, minimum_free_disk_gib=64,
        historical_405_startup_worker_seconds=sum(x['elapsed_seconds'] for x in earlier),
        four_roles_linear_startup_worker_seconds=measured,
        sensitivity_worker_seconds=[measured / 2, measured * 4],
        sum_wall_timeout_worker_hours=1620 * 900 / 3600,
        runtime_uncalibrated=True, finish_eta=None,
        available_memory_gib=psutil.virtual_memory().available / 2**30,
        available_disk_gib=shutil.disk_usage(out).free / 2**30,
        native_startup_only=True, posterior_sampling=False, gpu=False, new_cost_usd=0,
        caveat='Historical original startup, different initial state. Linear and sensitivity scenarios are planning figures; no memory or runtime guarantee. All capped/failed outcomes retained. Two simultaneous12GiB process address-space caps and24GiB cgroup constrain memory; no swap. Output allowance is an estimate, per-file cap enforced.')
    own = [Path(__file__), Path('scripts/baliphy_reference_initialization.py'),
           Path('scripts/run_baliphy_reference_preflight.py'), Path('scripts/readback_baliphy_reference_preflight.py'),
           Path('scripts/launch_baliphy_reference_preflight.py'), Path('scripts/check_baliphy_reference_initialization.py'),
           Path('scripts/check_baliphy_reference_preflight_grid.py'), Path('scripts/ancestral_chain_attempt.py'),
           Path('scripts/run_after_verified_dependencies_v2.py'), Path('scripts/record_project_runtime_checkpoint_v4.py'),
           Path('scripts/close_full_triad_sequence_stage.py'), Path('scripts/record_completed_process_handoffs_v2.py'),
           Path('scripts/reference_measurement_union_sources.py'), Path('scripts/run_ortholog_pair_guide_comparison.py'),
           a.software_validation, a.grid_validation, chain_path, proposals_path, closure_path, archive, source_manifest, binary, Path('/usr/bin/prlimit')]
    bindings.update({str(p): sha(p) for p in own + api_paths})
    for name, r in sources['files'].items(): bindings[str(source_manifest.parent / name)] = r['sha256']
    for group, row in models.items(): bindings[row['reference_program']] = row['reference_program_sha256']
    for path in [out / 'jobs.json', out / 'model_transformations.json']: bindings[str(path)] = sha(path)
    plan = dict(unit='fungal-baliphy-reference-preflight-20261003-v1.service', jobs=str(out / 'jobs.json'),
        inputs=str(out), output='results/ancestral/full-baliphy-reference-preflight-20261003-v1',
        completion='metadata/baliphy_reference_preflight_completed_20261003.json',
        completion_plan='metadata/baliphy_reference_preflight_completion_plan_20261003.json',
        launch_inventory='metadata/baliphy_reference_preflight_launches_20261003.json',
        dependencies=['metadata/baliphy_horizon_resources_closure_launch_20261003.json'],
        expected=dict(full_chains=1620, full_quartets=405, effective_inputs=135, original_configuration_aliases=324),
        software_validation=str(a.software_validation), resources=resources, pins=bindings, scope=SCOPE)
    with a.plan.open('x') as f: f.write(json.dumps(plan, indent=2) + '\n')
    print(json.dumps(dict(plan=str(a.plan), **resources)), flush=True)


if __name__ == '__main__': main()
