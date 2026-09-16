#!/usr/bin/env python3
"""Exercise the repair installation handoff in disposable synthetic directories."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(p, obj):
    p.write_text(json.dumps(obj))


def main():
    results = []
    for case in ['complete', 'wrong_tips', 'negative_branch', 'missing_length', 'occupied_target', 'missing_receipt']:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / 'source.fa'; source.write_text('>a\nAAAA\n>b\nCCCC\n>c\nDDDD\n')
            production = root / 'production.txt'
            production.write_text('retained existing data' if case == 'occupied_target' else '')
            previous = production.read_bytes()
            rebuild = root / 'rebuild'; rebuild.mkdir()
            trees = {'wrong_tips': '(a:1,b:2,d:3);', 'negative_branch': '(a:-1,b:2,c:3);',
                     'missing_length': '(a,b:2,c:3);'}
            (rebuild / 'tree.nwk').write_text(trees.get(case, '(a:1,b:2,c:3);'))
            (rebuild / 'fasttree.log').write_text('Synthetic fixture only\n')
            rebuild_script = 'scripts/rebuild_audited_family_tree.py'
            plan = root / 'plan.json'
            save(plan, dict(family='OGTEST', source=str(source), output=str(rebuild),
                            production_tree=str(production), pinned_files={
                                str(source): sha(source), str(production): sha(production),
                                rebuild_script: sha(ROOT / rebuild_script)}))
            save(rebuild / 'config.json', dict(plan_sha256=sha(plan), script_sha256=sha(ROOT / rebuild_script)))
            if case != 'missing_receipt':
                save(rebuild / 'receipt.json', dict(status='complete_isolated_audited_family_tree_rebuild',
                     family='OGTEST', tips=3, config_sha256=sha(rebuild / 'config.json'),
                     artifacts={n: sha(rebuild / n) for n in ['tree.nwk', 'fasttree.log']}))
            config = root / 'config.json'
            import os
            save(config, dict(producer={'pid': 2147483647, 'created': 0, 'command': []},
                 pins={str(plan): sha(plan)}, rebuild_plan=str(plan), output=str(root / 'out'),
                 cpu_affinity=[min(os.sched_getaffinity(0))], minimum_available_memory_bytes=0,
                 minimum_free_disk_bytes=0, python=str(ROOT / '.cache/envs/orthofinder/bin/python')))
            run = subprocess.run(['python', 'scripts/advance_repaired_family_installation.py',
                                  '--config', str(config)], cwd=ROOT, capture_output=True, text=True)
            state = json.loads((root / 'out/state.json').read_text())
            if case == 'complete':
                assert run.returncode == 0 and state['status'] == 'complete_verified_family_tree_installation'
                assert production.read_bytes() == (rebuild / 'tree.nwk').read_bytes()
            else:
                assert run.returncode != 0 and state['status'] == 'failed_requires_review'
                assert production.read_bytes() == previous
            results.append(dict(case=case, exit_code=run.returncode, status=state['status'],
                                production_preserved_on_failure=case != 'complete'))
    paths = ['scripts/advance_repaired_family_installation.py', 'scripts/validate_repaired_tree_native.py',
             'scripts/check_family_installation_gates.py']
    result = dict(status='passed_synthetic_family_installation_gates', cases=results,
                  scripts={p: sha(ROOT / p) for p in paths},
                  scope='Three-tip disposable fixtures; actual native parser and installer, including atomic replacement. Does not establish large-tree performance or biological validity.')
    out = ROOT / 'metadata/orthology_family_installation_gate_checks.json'
    out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
