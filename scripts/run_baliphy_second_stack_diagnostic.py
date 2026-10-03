#!/usr/bin/env python3
"""Preserve separate initial-state debugger attempts for the second failed input."""
import argparse
import copy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys

import psutil

from ancestral_chain_attempt import run_attempt, sha, write_json
from reference_measurement_union_sources import verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', type=int, choices=[4, 5], required=True)
    p.add_argument('--prepare', action='store_true')
    a = p.parse_args()
    root = Path(f'data/software_audits/baliphy-native-segfault-diagnostic-20261003-v{a.version}').resolve()
    cp = root / 'configuration.json'
    if a.prepare:
        root.mkdir(exist_ok=False)
        pp = Path('metadata/baliphy_joint_sampler_qualification_v3_plan_20261003.json')
        plan = json.loads(pp.read_text()); verify(plan['pins'])
        job = next(j for j in json.loads(Path(plan['jobs']).read_text())
                   if j['chain']['effective_input_group'].startswith('f9afda0e')
                   and j['chain']['prior_label'] == 'broad' and j['chain']['chain'] == 1)
        source = Path(plan['output']) / 'chains' / (job['chain']['chain_id'] + '.json')
        row = json.loads(source.read_text())
        assert row['exit_code'] == -11 and row['saved_alignments'] == 0
        command = copy.deepcopy(job['config']['command'])
        command[command.index('--iterations') + 1] = '0'
        command.insert(command.index('--'), '--core=0')
        if a.version == 5:
            command.insert(command.index('--'), '--stack=67108864')
        gdb = root / 'commands.gdb'
        gdb.write_text('set pagination off\nset confirm off\nset print thread-events off\nrun\nthread apply all bt 16\n' +
                       ('info registers rsp\ninfo proc mappings\n' if a.version == 4 else ''))
        config = dict(command=['/usr/bin/gdb', '--batch', '-x', str(gdb), '--args', *command],
                      timeout_seconds=900, pins=dict(job['config']['pins']))
        for fp in [gdb, Path('/usr/bin/gdb'), Path(__file__).resolve(), pp.resolve(), source.resolve()]:
            config['pins'][str(fp)] = sha(fp)
        write_json(cp, config)
        write_json(root / 'source_job.json', job)
        rp = Path(f'metadata/baliphy_native_segfault_diagnostic_resources_20261003_v{a.version}.json')
        assert not rp.exists()
        write_json(rp, dict(checked_utc=datetime.now(timezone.utc).isoformat(), source_role=row['chain_id'],
            source_role_checkpoint=str(source), source_role_checkpoint_sha256=sha(source),
            configuration=str(cp), configuration_sha256=sha(cp), cpus=2, memory_gib=64, swap_gib=0,
            native_address_space_gib=48, native_cpu_seconds=job['config']['timeout_seconds'],
            wall_seconds=900, service_wall_seconds=960, native_file_limit_gib=2, native_iterations=0,
            native_stack_bytes=64*2**20 if a.version == 5 else None, core_bytes=0, blas_threads=1,
            planning_output_gib=4, minimum_free_disk_gib=64, new_cost_usd=0, gpu=False,
            available_memory_gib=psutil.virtual_memory().available/2**30,
            available_disk_gib=shutil.disk_usage('.').free/2**30,
            scope='Isolated second-input representative initial-state causal diagnostic. Original model/input/tree/seed preserved; no posterior use, global default change or original restart. Gdb exit zero does not prove native success.'))
        print(str(cp), flush=True)
        return
    config = json.loads(cp.read_text()); verify(config['pins'])
    assert shutil.disk_usage(root).free >= 64*2**30
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    cg = Path('/sys/fs/cgroup') / group.lstrip('/')
    caps = {k: (cg/k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert caps == {'cpu.max':'200000 100000','memory.max':str(64*2**30),'memory.swap.max':'0'}
    assert all(os.environ.get(k) == '1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'])
    proc = psutil.Process()
    launch = dict(pid=proc.pid, created=proc.create_time(), cmdline=proc.cmdline(),
                  invocation_id=os.environ['INVOCATION_ID'], actual_cgroup_limits=caps)
    assert not (root/'controller.json').exists(); write_json(root/'controller.json', launch)
    print('original_controller', json.dumps(launch), flush=True)
    assert not (root/'attempt').exists(), 'Diagnostic attempts never restart'
    receipt = run_attempt(root/'attempt', config)
    print('debugger_receipt', str(receipt), 'exit', json.loads(receipt.read_text())['exit_code'], flush=True)


if __name__ == '__main__': main()
