#!/usr/bin/env python3
"""Copy audited full reconciliation inputs without enabling an incomplete restart."""
import argparse
import csv
import hashlib
import json
import shutil
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace
from orthofinder.utils import files


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    for name, expected in plan['pins'].items():
        if sha(Path(name)) != expected:
            raise ValueError(f'Pin mismatch: {name}')
    out = Path(plan['output']).resolve()
    if out.exists():
        raise FileExistsError(out)
    if shutil.disk_usage(out.parent).free < plan['resources']['minimum_free_disk_gib'] * 2**30:
        raise RuntimeError('Insufficient disk space')
    sources = [Path(p).resolve() for p in plan['source_working_directories']]
    universe = json.loads(Path(plan['universe_receipt']).read_text())
    inventory_receipt = json.loads(Path(plan['inventory_receipt']).read_text())
    inventory = Path(plan['inventory_receipt']).parent / 'families.tsv'
    if sha(inventory) != inventory_receipt['artifacts']['families.tsv']:
        raise ValueError('Inventory changed')
    with inventory.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    valid = [r for r in rows if r['tree_status'] == 'tip_branch_checks_passed']
    valid_names = {r['family'] for r in valid}
    pending = [r['family'] for r in rows if int(r['source_sequences']) >= 3 and r['family'] not in valid_names]
    if len(valid) != 32637 or pending != ['OG0000017']:
        raise ValueError('Unexpected tree disposition')
    species_file = sources[0] / 'SpeciesIDs.txt'
    species = [int(line.split(':', 1)[0]) for line in species_file.read_text().splitlines() if line.strip()]
    if len(species) != len(set(species)) or len(species) != 526:
        raise ValueError('Unexpected species IDs')
    if sha(species_file) != universe['species_ids_sha256']:
        raise ValueError('Species mapping changed')
    cluster = sources[0] / 'clusters_OrthoFinder.txt_id_pairs.txt'
    if sha(cluster) != universe['cluster_sha256']:
        raise ValueError('Cluster partition changed')
    inputs = [('SpeciesIDs.txt', species_file, sha(species_file)),
              ('SequenceIDs.txt', sources[0] / 'SequenceIDs.txt', sha(sources[0] / 'SequenceIDs.txt')),
              (cluster.name, cluster, sha(cluster))]
    for sid in species:
        candidates = [d / f'Species{sid}.fa' for d in sources if (d / f'Species{sid}.fa').is_file()]
        if len(candidates) != 1:
            raise ValueError(f'Ambiguous or missing species FASTA {sid}')
        inputs.append((candidates[0].name, candidates[0], sha(candidates[0])))
    for r in valid:
        name = f"Trees_ids/{r['family']}.txt"
        inputs.append((name, sources[0] / name, r['tree_sha256']))
    out.mkdir()
    start = time.time()
    reports = []
    guide_receipt = json.loads(Path(plan['guide_receipt']).read_text())
    for guide in guide_receipt['guides']:
        name = guide['guide']
        base = out / name
        source = base / 'Source'
        wd = source / 'WorkingDirectory'
        (wd / 'Trees_ids').mkdir(parents=True)
        manifest = base / 'copied_files.tsv'
        nbytes = 0
        with manifest.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['relative_path', 'source', 'bytes', 'sha256'])
            for relative, original, expected in inputs:
                dest = wd / relative
                if sha(original) != expected:
                    raise ValueError(f'Changed source: {original}')
                shutil.copyfile(original, dest)
                if dest.is_symlink() or dest.stat().st_ino == original.stat().st_ino or sha(dest) != expected or sha(original) != expected:
                    raise ValueError(f'Copy verification failed: {original}')
                size = dest.stat().st_size
                nbytes += size
                writer.writerow([str(dest.relative_to(base)), str(original), size, expected])
            for input_name, dest in [('species_tree_ids.nwk', wd / 'SpeciesTree_unrooted_ids.txt'),
                                     ('species_tree_taxa.nwk', base / 'species_tree_taxa.nwk')]:
                original = Path(plan['guide_directory']) / name / input_name
                expected = guide['artifacts'][input_name]
                if sha(original) != expected:
                    raise ValueError('Guide hash mismatch')
                shutil.copyfile(original, dest)
                if sha(dest) != expected:
                    raise ValueError('Guide copy mismatch')
                nbytes += dest.stat().st_size
                writer.writerow([str(dest.relative_to(base)), str(original.resolve()), dest.stat().st_size, expected])
        descriptor = (f'Isolated reconciliation inputs; pending validated OG0000017 installation.\n'
                      f'WorkingDirectory_Base: {wd}/\nFN_Orthogroups: {wd / cluster.name}\n'
                      f'WorkingDirectory_Trees: {wd}/\n')
        (source / 'Log.pending.txt').write_text(descriptor)
        # Probe the actual native locator in a disposable descriptor directory.
        # No executable Log.txt is left in the incomplete staged Source.
        with tempfile.TemporaryDirectory(prefix='native-locator-') as tmp:
            (Path(tmp) / 'Log.txt').write_text(descriptor)
            loc = files.PreviousFilesLocator_new(SimpleNamespace(qStartFromFasta=False, qStartFromBlast=False), tmp)
            if ([Path(p) for p in loc.wd_base_prev] != [wd]
                    or Path(loc.wd_trees) != wd
                    or Path(loc.clustersFilename_pairs) != wd / cluster.name):
                raise ValueError('Native locator escaped isolated inputs')
        reports.append(dict(guide=name, copied_files=len(inputs)+2, copied_bytes=nbytes,
                            taxa=len(species), validated_gene_trees=len(valid),
                            manifest_sha256=sha(manifest), pending_descriptor_sha256=sha(source / 'Log.pending.txt'),
                            native_locator_isolated=True, runnable_log_present=(source / 'Log.txt').exists()))
        print(name, reports[-1], flush=True)
    receipt = dict(status='staged_incomplete_full_reconciliation_inputs', guides=reports,
                   pending_families=pending, plan_sha256=sha(args.plan), script_sha256=sha(Path(__file__)),
                   elapsed_seconds=time.time()-start,
                   launch_ready=False,
                   remaining_gates=['Independently validated repaired OG0000017 copied into both variants',
                                    'Final complete input inventory and native restart resource plan',
                                    'Output validation beyond process exit status'],
                   interpretation='Independent physical copies for all 526 taxa and both conditional homogeneous guide alternatives. No live source files modified, no inference launched, and no runnable restart log published. Alignments/search intermediates are omitted for the previously tested --from-trees --no-fix-files path; runtime source-access coverage still requires full-run verification.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    main()
