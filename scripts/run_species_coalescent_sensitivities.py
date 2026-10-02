#!/usr/bin/env python3
"""Run all 30 supported ASTRAL-III candidate species-tree sensitivities.

Requires actual full input closure. Native annotated estimates remain
conditional on inferred gene trees, the MSC/locality model and support/taxon
settings. They do not establish roots, times, physical structural change or
biological causes of discordance. Full independent scoring is still required.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import dendropy
import psutil

from species_coalescent_sources import load
from reference_measurement_union_sources import bind, verify
from run_ortholog_pair_guide_comparison import sha


def write(path, record):
    with Path(path).open('x') as handle: handle.write(json.dumps(record, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    input_plan = json.loads(Path(plan['input_plan']).read_text())
    _, _, cohorts, _, bindings = load(input_plan, plan['input_plan'])
    bind(bindings, args.plan)
    for path, digest in plan['pins'].items(): bind(bindings, path, digest)
    completion = json.loads(Path(plan['input_completion']).read_text())
    assert completion['status'] == 'complete_verified_full_species_coalescent_inputs'
    assert completion['cases'] == 30 and completion['marker_states'] == 3750 and completion['exact_process_journals_checked'] == 2
    bind(bindings, plan['input_completion'])
    bind(bindings, completion['full_hash_archive'], completion['full_hash_archive_sha256'])
    archive = json.loads(Path(completion['full_hash_archive']).read_text())
    assert len(archive['services']) == 2
    for path, digest in archive['source_hashes'].items(): bind(bindings, path, digest)
    source = json.loads(Path(completion['producer_receipt']).read_text())
    assert len(source['summaries']) == source['cases'] == 30
    verify(bindings)
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'config.json', dict(plan_sha256=sha(args.plan), input_completion_sha256=sha(plan['input_completion']), source_hashes=bindings))
    bind(bindings, out / 'config.json')
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    runs = []
    for case in source['summaries']:
        folder = out / case['case']
        folder.mkdir()
        tree_path = folder / 'species.tree'
        command = [plan['java'], '-XX:ActiveProcessorCount=2', '-XX:+ExitOnOutOfMemoryError', '-Xmx28G',
                   '-jar', plan['jar'], '-i', str(Path(case['genes']).resolve()), '-o', str(tree_path.resolve()),
                   '-t', '2', '-s', str(plan['seed'])]
        with (folder / 'native.log').open('x') as handle:
            start = time.monotonic()
            process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT, env=env)
            actual = psutil.Process(process.pid)
            write(folder / 'native_launch.json', dict(pid=actual.pid, created=actual.create_time(), cmdline=actual.cmdline(), command=command, plan_sha256=sha(args.plan)))
            print('native_coalescent_case_started', case['case'], actual.pid, command, flush=True)
            code = process.wait()
        assert code == 0, (case['case'], code)
        assert tree_path.is_file()
        tree = dendropy.Tree.get(path=str(tree_path), schema='newick', rooting='force-unrooted', preserve_underscores=True)
        tips = [node.taxon.label for node in tree.leaf_node_iter()]
        assert len(tips) == len(set(tips)) == case['expected_taxa'] and set(tips) == cohorts[case['cohort']]['taxa']
        assert 'ASTRAL version 5.7.8' in (folder / 'native.log').read_text()
        record = dict(case=case, command=command, taxa=len(tips), returncode=code, elapsed_seconds=time.monotonic() - start,
                      artifacts={p.name: sha(p) for p in folder.iterdir() if p.is_file()})
        write(folder / 'receipt.json', dict(status='complete_native_coalescent_species_case_pending_full_independent_scoring', **record))
        bind(bindings, folder / 'receipt.json')
        for name, digest in record['artifacts'].items(): bind(bindings, folder / name, digest)
        runs.append(record)
        write(out / ('progress_' + str(len(runs)).zfill(2) + '.json'), dict(completed_cases=len(runs), total_cases=30, latest_case=case['case']))
        bind(bindings, out / ('progress_' + str(len(runs)).zfill(2) + '.json'))
        print('native_coalescent_case_completed', case['case'], len(runs), '/30', flush=True)
    assert len(runs) == 30
    verify(bindings)
    write(out / 'receipt.json', dict(status='complete_full_native_coalescent_species_sensitivities_pending_independent_scoring',
          plan_sha256=sha(args.plan), cases=30, alignments=2, cohorts=5, support_policies=3, markers_per_case=125,
          input_completion_sha256=sha(plan['input_completion']), runs=runs, source_hashes=bindings,
          scientific_eligibility=False, scope=plan['scope']))


if __name__ == '__main__':
    main()
