#!/usr/bin/env python3
"""Qualify real signal-stack capture and normal exit on two private targets."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

from ancestral_chain_attempt import sha
from run_scalar_native_signal_debugger_v1 import run
from reference_measurement_union_sources import bind, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    assert not a.receipt.exists()
    a.output.mkdir(exist_ok=False)
    pins = {str(path): sha(path) for path in [Path(__file__),
        Path('scripts/run_scalar_native_signal_debugger_v1.py'), Path(sys.executable),
        Path('/usr/bin/gdb'), Path('/usr/bin/prlimit')]}
    outcomes = {}
    for case, code in [('segfault', 'import ctypes;ctypes.string_at(0)'), ('normal', 'print("native-normal-exit")')]:
        target = a.output/(case+'.py')
        target.write_text(code+'\n')
        bind(pins, target)
        plan = dict(output=str(a.output/(case+'-debugger')), native_command=[sys.executable, str(target.resolve())],
            native_caps=dict(address_space_bytes=2**30, cpu_seconds=20, per_file_bytes=16*2**20),
            debugger_wall_seconds=30, pins=dict(pins), artificial_control=True,
            scope='Private artificial actual ctypes null-pointer native signal or normal exit; no biological run or repair.')
        pp = a.output/(case+'.plan.json')
        with pp.open('x') as handle:json.dump(plan, handle, indent=2)
        rp = a.output/(case+'.receipt.json')
        result = run(pp, rp)
        assert result['native_stop']['live_stop'] == (case == 'segfault')
        if case == 'segfault':
            assert result['native_stop']['signal'] == 'SIGSEGV'
            assert '#0 ' in Path(plan['output'], 'gdb_stdout.log').read_text()
        else:
            assert 'No registers.' not in Path(plan['output'], 'gdb_stderr.log').read_text()
        outcomes[case] = dict(status=result['status'], signal=result['native_stop']['signal'], receipt=str(rp), sha256=sha(rp))
    for path in a.output.rglob('*'):
        if path.is_file():bind(pins, path)
    verify(pins)
    proof = dict(status='passed_actual_native_signal_and_normal_exit_debugger_controls',
        checked_utc=datetime.now(timezone.utc).isoformat(), outcomes=outcomes, source_hashes=pins,
        artificial_controls=2, original_attempts_restarted=False, scientific_eligibility=False, gpu=False)
    with a.receipt.open('x') as handle:json.dump(proof, handle, indent=2);handle.write('\n')
    print(json.dumps({k: v for k, v in proof.items() if k != 'source_hashes'}, indent=2))


if __name__ == '__main__':main()
