#!/usr/bin/env python3
"""Capture a separate native signal stop with original command and process caps."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

from ancestral_chain_attempt import sha
from reference_measurement_union_sources import bind, verify


def run(plan_path, receipt_path):
    plan_path, receipt_path = Path(plan_path), Path(receipt_path)
    assert not receipt_path.exists()
    plan = json.loads(plan_path.read_text())
    pins = dict(plan['pins'])
    verify(pins)
    root = Path(plan['output']).resolve()
    root.mkdir(parents=True, exist_ok=False)
    command = plan['native_command']
    assert command and all(isinstance(v, str) and '\n' not in v for v in command)
    caps = plan['native_caps']
    script = root/'gdb_commands.txt'
    script.write_text('set pagination off\nset confirm off\n'
        'set disable-randomization off\nhandle SIGSEGV stop print nopass\n'
        'handle SIGABRT stop print nopass\nset exec-wrapper /usr/bin/prlimit '
        '--as='+str(caps['address_space_bytes'])+' --cpu='+str(caps['cpu_seconds'])+
        ' --fsize='+str(caps['per_file_bytes'])+' --\npython\n'
        'import gdb,json,os\nproject_signal=None\n'
        'def project_stopped(event):\n'
        '    global project_signal\n'
        '    if isinstance(event,gdb.SignalEvent):project_signal=event.stop_signal\n'
        'gdb.events.stop.connect(project_stopped)\nend\nrun\npython\n'
        'pid=gdb.selected_inferior().pid\n'
        'record={"pid":pid,"signal":project_signal,"live_stop":pid>0}\n'
        'if pid>0:\n'
        '    record["pc"]=str(gdb.parse_and_eval("$pc"))\n'
        '    record["native_cmdline"]=open("/proc/"+str(pid)+"/cmdline","rb").read().decode().rstrip("\\0").split("\\0")\n'
        '    for name in ["status","limits","maps"]:\n'
        '        text=open("/proc/"+str(pid)+"/"+name).read()\n'
        '        open("native-stop-"+name+".txt","x").write(text)\n'
        'print("PROJECT_NATIVE_STOP "+json.dumps(record,sort_keys=True))\n'
        'if pid>0:\n'
        '    for cmd in ["p $_siginfo","thread apply all bt 32","info registers","x/16i $pc-32","info sharedlibrary"]:\n'
        '        try:gdb.execute(cmd)\n'
        '        except gdb.error as error:print("PROJECT_INSPECTION_ERROR "+str(error))\n'
        'else:print("PROJECT_NORMAL_EXIT_INSPECTION_SKIPPED")\nend\n')
    debugger = ['/usr/bin/gdb', '-nx', '--batch', '-x', str(script), '--args', *command]
    stdout, stderr = root/'gdb_stdout.log', root/'gdb_stderr.log'
    with stdout.open('x') as out, stderr.open('x') as err:
        result = subprocess.run(debugger, cwd=root, stdout=out, stderr=err,
                                timeout=plan['debugger_wall_seconds'])
    text = stdout.read_text()
    records = [json.loads(line[len('PROJECT_NATIVE_STOP '):]) for line in text.splitlines()
               if line.startswith('PROJECT_NATIVE_STOP ')]
    assert result.returncode == 0 and len(records) == 1
    record = records[0]
    if record['live_stop']:
        assert record['signal'] in ['SIGSEGV', 'SIGABRT'] and record['native_cmdline'] == command
        limits = (root/'native-stop-limits.txt').read_text()
        line = next(line for line in limits.splitlines() if line.startswith('Max address space'))
        assert str(caps['address_space_bytes']) in line
    else:
        assert record['signal'] is None and 'PROJECT_NORMAL_EXIT_INSPECTION_SKIPPED' in text
    artifacts = {str(path): sha(path) for path in root.rglob('*') if path.is_file()}
    for path in [Path(__file__), plan_path]:
        bind(pins, path)
    verify(pins)
    proof = dict(status='captured_separate_native_signal_stop' if record['live_stop']
                 else 'separate_native_debugger_exit_without_signal_stop',
        checked_utc=datetime.now(timezone.utc).isoformat(), plan_sha256=sha(plan_path),
        native_stop=record, debugger_exit_code=result.returncode, native_caps=caps,
        debugger_command=debugger, artifacts=artifacts, source_hashes=pins,
        original_attempts_restarted=False, source_or_prior_changed=False,
        scientific_eligibility=False, posterior_qualified=False, gpu=False,
        artificial_control=plan['artificial_control'], scope=plan['scope'])
    with receipt_path.open('x') as handle:
        json.dump(proof, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in proof.items() if k not in ['artifacts', 'source_hashes']}, indent=2))
    return proof


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    run(args.plan, args.receipt)
