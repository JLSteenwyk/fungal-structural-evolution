#!/usr/bin/env python3
"""Check actual completed/failed original journals and reject false handoffs."""
import argparse
import copy
import json
from pathlib import Path
import subprocess

from run_after_verified_dependencies_v2 import completed,inspect_journal
from run_ortholog_pair_guide_comparison import sha


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    success_path='metadata/species_coalescent_native_launch_20261002.json'
    failure_path='metadata/species_coalescent_quartet_launch_20261002.json'
    good=json.loads(Path(success_path).read_text());bad=json.loads(Path(failure_path).read_text())
    actual=completed(good)
    # Collected failed unit can report default success fields; its original journal must reject it.
    try:completed(bad)
    except (RuntimeError,AssertionError):failed_rejected=True
    else:raise AssertionError('Original failed dependency accepted')
    raw=subprocess.check_output(['journalctl','--user','-u',good['unit'],'-o','json','--no-pager'],text=True)
    rows=[json.loads(line) for line in raw.splitlines()]
    invocation=next(r['_SYSTEMD_INVOCATION_ID'] for r in rows if r.get('_PID')==str(good['pid']))
    exact=[r for r in rows if r.get('_PID')==str(good['pid']) and r.get('_CMDLINE')==' '.join(good['cmdline'])]
    resources=[r for r in rows if r.get('USER_INVOCATION_ID')==invocation and r.get('CPU_USAGE_NSEC')]
    assert resources and exact
    cases=[('no_original_messages',resources),('no_original_resources',exact),
           ('foreign_resource_invocation',exact+[dict(r,USER_INVOCATION_ID='foreign') for r in resources]),
           ('changed_original_command',[dict(r,_CMDLINE='changed') for r in exact]+resources),
           ('failed_original_invocation',exact+resources+[dict(USER_INVOCATION_ID=invocation,MESSAGE='Failed with result exit-code')])]
    rejected=0
    for name,evidence in cases:
        try:inspect_journal(good,evidence)
        except (RuntimeError,AssertionError):rejected+=1
        else:raise AssertionError('False journal handoff accepted '+name)
    inspect_journal(good,exact+resources+[dict(USER_INVOCATION_ID=invocation,MESSAGE=None)])
    result=dict(status='passed_original_dependency_completion_journal_contracts',
        actual_native_completion_accepted=True,actual_failed_auditor_rejected=failed_rejected,
        altered_journals_rejected=rejected,null_message_accepted=True,
        original_native_journal_sha256=__import__('hashlib').sha256(raw.encode()).hexdigest(),
        source_hashes={p:sha(p) for p in [str(Path(__file__)),'scripts/run_after_verified_dependencies_v2.py',
            'scripts/run_after_verified_dependencies.py',success_path,failure_path]},scientific_eligibility=False,
        scope='Actual original30-case native completion and actual failed full auditor invocation checked. No-message/no-resource/foreign-invocation/changed-command/failure and null-message journal contracts; collected systemd default success/0 alone never proves completion. Software handoff check only; production numerical and biological acceptance remain separate.')
    with (args.output/'receipt.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if isinstance(v,(str,int,bool))}),flush=True)


if __name__=='__main__':main()
