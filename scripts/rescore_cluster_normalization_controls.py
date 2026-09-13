#!/usr/bin/env python3
"""Rescore unchanged saved alignments with matched original/patched Foldseek builds."""
import argparse
import json
import subprocess
from pathlib import Path
from audit_busco_gene_copies import sha
from assess_pae_sensitivity import checked_receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['build', 'exact', 'output']:
        p.add_argument('--' + name, required=True, type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('Use new immutable output')
    checked_receipt(a.build)
    config = json.loads((a.exact / 'config.json').read_text())
    completion = json.loads((a.exact / 'completion.json').read_text())
    if completion['config_sha256'] != sha(a.exact / 'config.json') or completion['alignment_table_sha256'] != sha(a.exact / 'alignments.tsv'):
        raise ValueError('Exact source receipt mismatch')
    command = config['command']
    if command[command.index('--exact-tmscore') + 1] != '1':
        raise ValueError('Expected exact score protocol')
    a.output.mkdir(parents=True)
    # Input database components, including the saved alignment DB, must remain fixed.
    pins = {}
    for prefix in set(command[2:5]):
        path = Path(prefix)
        candidates = [f for f in path.parent.glob(path.name + '*') if f.is_file() and (f.name == path.name or f.name.startswith(path.name + '.') or f.name.startswith(path.name + '_'))]
        if not candidates:
            raise ValueError('Missing source database')
        for f in candidates:
            pins[str(f)] = sha(f)
    plan = {'build_receipt_sha256': sha(a.build / 'receipt.json'),
            'exact_completion_sha256': sha(a.exact / 'completion.json'),
            'script_sha256': sha(Path(__file__)), 'database_files': pins,
            'cpu_threads': 8, 'memory_gb_planning': 32, 'disk_gb_planning': 5,
            'wall_hours_planning': [.5, 12],
            'estimate_basis': 'Prior identical-size exact conversion required about 17.5 minutes; two sequential control conversions, allowing toolchain variation.',
            'cost': 'Existing authorized host; no new charges', 'commands': []}
    for variant in ['unmodified', 'inclusive-span']:
        cmd = list(command)
        cmd[0] = str((a.build / ('foldseek-' + variant)).resolve())
        cmd[5] = str(a.output / (variant + '.tsv'))
        plan['commands'].append(cmd)
    (a.output / 'config.json').write_text(json.dumps(plan, indent=2) + '\n')
    for variant, cmd in zip(['unmodified', 'inclusive-span'], plan['commands']):
        with (a.output / (variant + '.log')).open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        (a.output / (variant + '-completion.json')).write_text(json.dumps({'status': 'conversion_complete_pending_readback', 'command': cmd, 'output_sha256': sha(Path(cmd[5]))}, indent=2) + '\n')
    for path, digest in pins.items():
        if sha(Path(path)) != digest:
            raise ValueError('Database changed during conversion')
    result = {'status': 'complete_normalization_control_conversions_pending_readback',
              'config_sha256': sha(a.output / 'config.json'),
              'interpretation': 'Two matched-toolchain conversions on unchanged alignments. No corrected membership or accuracy conclusion before full numerical readback.',
              'artifacts': {f.name: sha(f) for f in a.output.iterdir() if f.suffix == '.tsv'}}
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
