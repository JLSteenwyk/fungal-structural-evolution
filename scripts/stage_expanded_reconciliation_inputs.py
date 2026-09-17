#!/usr/bin/env python3
"""Assemble isolated expanded inputs, preserving explicit missing-tree gates."""
import argparse
import csv
import hashlib
import json
import shutil
import time
from pathlib import Path


def sha(p):
    with p.open('rb') as h:
        return hashlib.file_digest(h, 'sha256').hexdigest()


def read(p):
    return json.loads(p.read_text())


def rows(p):
    with p.open() as h:
        return list(csv.DictReader(h, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    a = ap.parse_args()
    plan = read(a.plan)
    for name, digest in plan['pins'].items():
        if sha(Path(name)) != digest:
            raise ValueError('Changed pin: ' + name)
    output = Path(plan['output']).resolve()
    if output.exists():
        raise FileExistsError(output)
    if shutil.disk_usage(output.parent).free < plan['resources']['minimum_free_disk_gib']*2**30:
        raise RuntimeError('Insufficient free disk')
    staged, merged, retained, new = [Path(plan[k]) for k in ['staged', 'merged', 'retained', 'new_trees']]
    sr, mr, rr, nr = [read(p) for p in [staged/'receipt.json', merged/'receipt.json', retained/'receipt.json', new/'readback.json']]
    if nr['status'] != 'passed_complete_native_new_family_tree_readback':
        raise ValueError('New trees have not passed native readback')
    if sha(retained/'retained_tree_catalog.tsv') != rr['catalog_sha256']:
        raise ValueError('Retained catalog changed')
    catalog = {r['membership_sha256']: r for r in rows(retained/'retained_tree_catalog.tsv')}
    for name, digest in nr['artifacts'].items():
        if sha(new/name) != digest:
            raise ValueError('Native readback artifact changed')
    new_receipts = read(new/'native_readback_family_receipts.json')
    output.mkdir()
    started = time.time()
    reports = []
    for guide in mr['guides']:
        name = guide['guide']
        for relative, digest in guide['artifacts'].items():
            if sha(merged/relative) != digest:
                raise ValueError('Merged artifact changed')
        prior = next(g for g in sr['guides'] if g['guide']==name)
        manifest_path = staged/name/'copied_files.tsv'
        if sha(manifest_path) != prior['manifest_sha256']:
            raise ValueError('Staged manifest changed')
        base = output/name
        wd = base/'Source/WorkingDirectory'
        (wd/'Trees_ids').mkdir(parents=True)
        manifest = base/'copied_files.tsv'
        total_bytes = files = trees = 0
        pending = []
        with manifest.open('w') as h:
            writer = csv.writer(h, delimiter='\t', lineterminator='\n')
            writer.writerow(['relative_path','source','bytes','sha256'])

            def copy(source, relative, digest):
                nonlocal total_bytes, files
                dest = base/relative
                if dest.exists() or sha(source) != digest:
                    raise ValueError('Existing destination or changed source: ' + str(source))
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, dest)
                if dest.is_symlink() or dest.stat().st_ino == source.stat().st_ino or sha(dest) != digest or sha(source) != digest:
                    raise ValueError('Physical copy verification failed')
                size = dest.stat().st_size
                total_bytes += size
                files += 1
                writer.writerow([relative,str(source.resolve()),size,digest])

            for row in rows(manifest_path):
                relative = row['relative_path']
                if '/Trees_ids/' in relative or relative.endswith('clusters_OrthoFinder.txt_id_pairs.txt'):
                    continue
                copy(staged/name/relative, relative, row['sha256'])
            cluster = merged/name/'clusters_id_pairs.txt'
            copy(cluster, 'Source/WorkingDirectory/clusters_OrthoFinder.txt_id_pairs.txt', sha(cluster))
            crosswalk = merged/name/'family_sources.tsv'
            copy(crosswalk, 'family_sources.tsv', sha(crosswalk))
            families = rows(crosswalk)
            if len(families) != guide['families'] or sum(int(r['proteins']) for r in families) != 5815847:
                raise ValueError('Incomplete expanded partition')
            for row in families:
                if int(row['proteins']) < 3:
                    continue
                key = row['membership_sha256']
                if row['source_type'] == 'retained':
                    old = catalog[key]
                    if old[name+'_family'] != row['new_family'] or old['original_family'] != row['original_tree_candidate']:
                        raise ValueError('Retained destination mismatch')
                    if old['status'] == 'pending_missing_staged_tree':
                        pending.append(dict(family=row['new_family'],original_family=old['original_family'],membership_sha256=key))
                        continue
                    if old['status'] != 'validated_exact_membership':
                        raise ValueError('Unvalidated retained tree')
                    source, digest = Path(old['tree_path']), old['tree_sha256']
                elif row['source_type'] == 'discovery':
                    folder = new/'families'/key
                    if sha(folder/'receipt.json') != new_receipts[key]:
                        raise ValueError('New family receipt differs from native readback')
                    receipt = read(folder/'receipt.json')
                    if receipt['membership_sha256'] != key or receipt['dimensions']['proteins'] != int(row['proteins']):
                        raise ValueError('New family membership mismatch')
                    source, digest = folder/'tree.nwk', receipt['artifacts']['tree.nwk']
                else:
                    raise ValueError('Unexpected source type')
                copy(source, 'Source/WorkingDirectory/Trees_ids/'+row['new_family']+'.txt', digest)
                trees += 1
        if [p['original_family'] for p in pending] != ['OG0000017']:
            raise ValueError('Unexpected pending families')
        descriptor = (f'Expanded inputs pending repaired family and independent input readback.\n'
                      f'WorkingDirectory_Base: {wd}/\n'
                      f'FN_Orthogroups: {wd}/clusters_OrthoFinder.txt_id_pairs.txt\n'
                      f'WorkingDirectory_Trees: {wd}/\n')
        (base/'Source/Log.pending.txt').write_text(descriptor)
        report = dict(guide=name,families=len(families),proteins=5815847,copied_trees=trees,
                      pending_families=pending,copied_files=files,copied_bytes=total_bytes,
                      manifest_sha256=sha(manifest),pending_descriptor_sha256=sha(base/'Source/Log.pending.txt'))
        reports.append(report)
        print(json.dumps(report),flush=True)
    for name, digest in plan['pins'].items():
        if sha(Path(name)) != digest:
            raise ValueError('Pin changed during staging')
    result = dict(status='staged_expanded_inputs_requires_readback_and_repaired_tree',
                  launch_ready=False,guides=reports,plan_sha256=sha(a.plan),
                  script_sha256=sha(Path(__file__)),elapsed_seconds=time.time()-started,
                  scope='Complete expanded partitions and all available validated trees copied using exact family mappings; singleton and pair families retained. No runnable Log.txt, no reconciliation launched. Both guides remain conditional homogeneous alternatives.')
    (output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
