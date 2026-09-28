#!/usr/bin/env python3
"""Diagnose each fully checked quartet; retain unresolved terminal groups."""
import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time

import psutil
from ancestral_chain_attempt import sha, write_json


def quartet_manifest(jobs, root, iterations, mapping_hash):
    available = [list((root / job['chain']['chain_id']).glob('attempt-*-sample-audit.json')) for job in jobs]
    if any(not paths for paths in available):
        return None
    assert all(len(paths) == 1 for paths in available)
    chains, evidence = [], {}
    for job, paths in zip(jobs, available):
        chain, config = job['chain'], job['config']
        folder = root / chain['chain_id']
        ap = paths[0]
        audit = json.loads(ap.read_text())
        assert audit['status'] == 'all_saved_alignments_and_candidate_nodes_checked'
        assert audit['iterations'] == iterations and audit['mapping_sha256'] == mapping_hash
        rp = Path(audit['attempt_receipt'])
        assert rp.resolve().parent.parent == folder.resolve()
        assert sha(rp) == audit['attempt_receipt_sha256']
        receipt = json.loads(rp.read_text())
        assert receipt['status'] == 'exited_zero_pending_scientific_validation' and receipt['exit_code'] == 0
        digest = hashlib.sha256(json.dumps(config, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        assert receipt['configuration_sha256'] == digest
        assert json.loads((folder / 'configuration.json').read_text()) == config
        assert json.loads((rp.parent / 'command.json').read_text()) == config['command']
        for name, expected in receipt['artifacts'].items():
            assert sha(rp.parent / name) == expected
        assert sha(audit['scalar_log']) == audit['scalar_log_sha256']
        assert config['seed'] == chain['seed']
        chains.append(dict(chain_id=chain['chain_id'], seed=chain['seed'],
            model_input_identity=config['model_input_identity'],
            log=audit['scalar_log'], log_sha256=audit['scalar_log_sha256']))
        evidence[str(ap)] = sha(ap)
        evidence[str(rp)] = sha(rp)
    assert len(chains) == len({c['seed'] for c in chains}) == 4
    assert len({c['model_input_identity'] for c in chains}) == 1
    headers = []
    for chain in chains:
        with open(chain['log']) as handle:
            headers.append(next(csv.reader(handle, delimiter='\t')))
    assert all(h == headers[0] for h in headers)
    fixed = {'iter', 'scale', 'scale1', 'scale*|T|', 'scale1*|T|', '|T|'}
    assert fixed <= set(headers[0])
    variables = [name for name in headers[0] if name not in fixed]
    assert {'prior', 'likelihood', 'posterior', 'ASRV.Gamma:alpha', 'RS07:rate', 'RS07:meanLength'} <= set(variables)
    assert len([name for name in variables if name.startswith('F:pi[')]) == 20
    return dict(chains=chains, variables=variables, expected_iterations=list(range(iterations + 1))), evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    plan_hash = sha(args.plan)
    def verify():
        assert sha(args.plan) == plan_hash
        for path, digest in plan['pins'].items():
            assert sha(path) == digest, path
    verify()
    source = json.loads(Path(plan['producer_plan']).read_text())
    launch = json.loads(Path(plan['producer_launch']).read_text())
    assert launch['plan_sha256'] == sha(plan['producer_plan'])
    jobs = json.loads(Path(source['jobs']).read_text())
    groups = defaultdict(list)
    for job in jobs:
        groups[job['config']['model_input_identity']].append(job)
    assert len(groups) == 405 and all(len(jobs) == 4 for jobs in groups.values())
    out = Path(plan['output'])
    out.mkdir(parents=True, exist_ok=False)
    root = Path(source['output'])
    complete = {}
    while True:
        verify()
        state = dict(line.split('=', 1) for line in subprocess.check_output([
            'systemctl', '--user', 'show', launch['unit'], '-p', 'MainPID', '-p', 'ActiveState',
            '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        terminal = state['ActiveState'] in ['inactive', 'failed']
        if not terminal:
            assert state['ActiveState'] == 'active' and int(state['MainPID']) == launch['pid']
            try:
                process = psutil.Process(launch['pid'])
                assert process.create_time() == launch['created'] and process.cmdline() == launch['cmdline']
            except psutil.NoSuchProcess:
                time.sleep(1)
                continue
        for group, group_jobs in sorted(groups.items()):
            if group in complete:
                continue
            ready = quartet_manifest(group_jobs, root, source['iterations'], sha(source['mapping']))
            if ready is None:
                continue
            manifest, evidence = ready
            folder = out / group
            folder.mkdir()
            outputs = {}
            for fraction in source['diagnostics']['burn_in_fractions']:
                cutoff = int(source['iterations'] * fraction)
                manifest['discard_through_iteration'] = cutoff
                mp = folder / ('manifest-%d.json' % cutoff)
                write_json(mp, manifest)
                destination = folder / ('discard-%d' % cutoff)
                subprocess.run([plan['python'], plan['diagnostic_script'], '--manifest', str(mp),
                                '--output', str(destination)], check=True, stdout=subprocess.DEVNULL)
                dp = destination / 'diagnostics.json'
                report = json.loads(dp.read_text())
                assert report['status'] == 'scalar_diagnostics_complete_not_posterior_qualification'
                assert set(report['variables']) == set(manifest['variables'])
                for path, digest in report['pins'].items():
                    assert sha(path) == digest
                outputs[str(cutoff)] = dict(path=str(dp), sha256=sha(dp),
                    counts=dict(Counter(r['status'] for r in report['variables'].values())))
            result = dict(status='both_burnin_scalar_screens_complete', evidence=evidence, outputs=outputs,
                          scope='Scalar screens only; no ancestral/alignment convergence qualification.')
            write_json(folder / 'receipt.json', result)
            complete[group] = dict(status=result['status'], receipt=str(folder / 'receipt.json'),
                                   receipt_sha256=sha(folder / 'receipt.json'))
            print('Diagnosed quartets', len(complete), '/405', group, flush=True)
        if terminal:
            break
        time.sleep(30)
    dispositions = {group: complete.get(group, dict(status='unresolved_missing_checked_quartet',
                    chain_ids=[j['chain']['chain_id'] for j in group_jobs])) for group, group_jobs in groups.items()}
    write_json(out / 'receipt.json', dict(status='terminal_group_accounting_requires_scientific_review',
        producer_terminal_state=state, plan_sha256=plan_hash, groups=dispositions,
        complete_quartets=len(complete), unresolved_quartets=405-len(complete),
        scope='Retains incomplete groups and producer failures. No global posterior qualification.'))


if __name__ == '__main__':
    main()
