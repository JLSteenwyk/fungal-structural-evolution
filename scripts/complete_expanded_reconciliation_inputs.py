#!/usr/bin/env python3
"""Create fresh expanded inputs including the independently repaired final tree."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from orthofinder.tools import tree


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    plan = read(args.plan)
    plan_sha = sha(args.plan)

    def verify():
        if sha(args.plan) != plan_sha:
            raise ValueError('Plan changed')
        for name, expected in plan['pins'].items():
            if sha(Path(name)) != expected:
                raise ValueError('Changed pin: ' + name)

    verify()
    source, output = [Path(plan[k]).resolve() for k in ['source', 'output']]
    if output.exists():
        raise FileExistsError(output)
    if shutil.disk_usage(output.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise ValueError('Insufficient storage')
    previous = read(source / 'receipt.json')
    previous_audit = read(source / 'readback.json')
    if (previous_audit['status'] != 'passed_complete_available_expanded_input_readback'
            or previous_audit['input_receipt_sha256'] != sha(source / 'receipt.json')):
        raise ValueError('Previous inputs lack bound independent readback')
    installation = read(Path(plan['installation']))
    native = read(Path(plan['native_readback']))
    repaired = Path(installation['production_tree'])
    if (installation['status'] != 'complete_verified_family_tree_installation'
            or native['status'] != 'passed_native_repaired_tree_readback'
            or sha(Path(plan['native_readback'])) != installation['native_readback_sha256']
            or sha(repaired) != installation['installed_tree_sha256']
            or sha(repaired) != native['tree_sha256']):
        raise ValueError('Unverified repaired tree')
    tips = tree.Tree(repaired.read_text().strip(), format=0).get_leaf_names()
    if len(tips) != len(set(tips)) or len(tips) != native['tips']:
        raise ValueError('Repaired tree tip universe differs')
    membership = hashlib.sha256(('\n'.join(sorted(tips)) + '\n').encode()).hexdigest()
    output.mkdir()
    started = time.time()
    reports = []
    for guide in previous['guides']:
        name = guide['guide']
        old_base, base = source / name, output / name
        pending = guide['pending_families']
        if (len(pending) != 1 or pending[0]['original_family'] != 'OG0000017'
                or pending[0]['membership_sha256'] != membership):
            raise ValueError('Repair does not match the pending family')
        old_manifest = old_base / 'copied_files.tsv'
        if sha(old_manifest) != guide['manifest_sha256']:
            raise ValueError('Old manifest changed')
        copied = rows(old_manifest)
        crosswalk = rows(old_base / 'family_sources.tsv')
        mapped = [r for r in crosswalk if r['new_family'] == pending[0]['family']]
        if (len(mapped) != 1 or mapped[0]['membership_sha256'] != membership
                or mapped[0]['source_type'] != 'retained'
                or mapped[0]['original_tree_candidate'] != 'OG0000017'
                or int(mapped[0]['proteins']) != len(tips)):
            raise ValueError('Repaired tree crosswalk differs')
        paths = [r['relative_path'] for r in copied]
        repair_relative = 'Source/WorkingDirectory/Trees_ids/' + pending[0]['family'] + '.txt'
        if len(set(paths)) != len(paths) or repair_relative in paths:
            raise ValueError('Duplicate or already installed tree path')
        base.mkdir()
        size = 0
        with (base / 'copied_files.tsv').open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['relative_path', 'source', 'bytes', 'sha256'])

            def copy(src, relative, digest):
                nonlocal size
                dest = base / relative
                if dest.exists() or sha(src) != digest:
                    raise ValueError('Existing target or changed source')
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, dest)
                if (dest.is_symlink() or dest.stat().st_ino == src.stat().st_ino
                        or sha(dest) != digest or sha(src) != digest):
                    raise ValueError('Physical copy mismatch')
                size += dest.stat().st_size
                writer.writerow([relative, str(src.resolve()), dest.stat().st_size, digest])

            for row in copied:
                copy(old_base / row['relative_path'], row['relative_path'], row['sha256'])
            copy(repaired, repair_relative, native['tree_sha256'])
        wd = base / 'Source/WorkingDirectory'
        descriptor = base / 'Source/Log.pending.txt'
        descriptor.write_text('Complete tree inputs; reconciliation launch remains separately gated.\n'
                              f'WorkingDirectory_Base: {wd}/\n'
                              f'FN_Orthogroups: {wd}/clusters_OrthoFinder.txt_id_pairs.txt\n'
                              f'WorkingDirectory_Trees: {wd}/\n')
        reports.append(dict(guide=name, families=guide['families'], proteins=guide['proteins'],
                            copied_trees=guide['copied_trees'] + 1, pending_families=[],
                            copied_files=len(copied) + 1, copied_bytes=size,
                            manifest_sha256=sha(base / 'copied_files.tsv'),
                            pending_descriptor_sha256=sha(descriptor)))
        print(json.dumps(reports[-1]), flush=True)
    verify()
    receipt = dict(status='staged_complete_expanded_inputs_requires_independent_readback',
                   guides=reports, launch_ready=False, plan_sha256=plan_sha,
                   previous_receipt_sha256=sha(source / 'receipt.json'),
                   repaired_tree_sha256=native['tree_sha256'],
                   script_sha256=sha(Path(__file__)), elapsed_seconds=time.time() - started,
                   scope='Fresh physical copies preserve prior snapshots. All candidate tree files staged; reconciliation and species-tree validation remain separate.')
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    with (output / 'readback.log').open('w') as log:
        subprocess.run([sys.executable, 'scripts/readback_expanded_reconciliation_inputs.py',
                        '--inputs', str(output), '--output', str(output / 'readback.json')],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    verify()
    audit = read(output / 'readback.json')
    if (audit['status'] != 'passed_complete_expanded_input_readback'
            or any(g['pending'] for g in audit['guides'])
            or audit['input_receipt_sha256'] != sha(output / 'receipt.json')):
        raise ValueError('Independent final readback failed')
    print('Complete input staging and native readback passed', flush=True)


if __name__ == '__main__':
    main()
