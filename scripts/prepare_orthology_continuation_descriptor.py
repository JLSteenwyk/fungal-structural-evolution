#!/usr/bin/env python3
"""Prepare and natively read back an explicit descriptor for interrupted orthology."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from orthofinder.utils import files


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def locate(folder):
    options = SimpleNamespace(qStartFromFasta=False, qStartFromBlast=False)
    locator = files.PreviousFilesLocator_new(options, str(folder.resolve()))
    return dict(working_directory_chain=locator.wd_base_prev,
                clusters=locator.clustersFilename_pairs, tree_working_directory=locator.wd_trees,
                rooted_species_tree=locator.speciesTreeRootedIDsFN,
                unrooted_species_tree=locator.speciesTreeUNRootedIDsFN)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--universe-receipt', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    source = a.source.resolve()
    wd = source / 'WorkingDirectory'
    log = source / 'Log.txt'
    clusters = wd / 'clusters_OrthoFinder.txt_id_pairs.txt'
    ids = wd / 'SpeciesIDs.txt'
    universe = json.loads(a.universe_receipt.read_text())
    if (universe['status'] != 'passed_complete_cluster_to_family_identity_readback'
            or sha(clusters) != universe['cluster_sha256']
            or sha(ids) != universe['species_ids_sha256']):
        raise ValueError('Cluster identity audit mismatch')
    log_hash = sha(log)
    before = locate(source)
    original = log.read_text()
    if before['clusters'] is not None or any(line.startswith('WorkingDirectory_Trees: ') for line in original.splitlines()):
        raise ValueError('Source metadata differs from the reviewed interrupted state')
    if not (wd / 'Trees_ids').is_dir():
        raise ValueError('No source gene-tree directory')
    missing_tree = wd / 'Trees_ids/OG0000017.txt'
    if not missing_tree.exists():
        raise ValueError('Expected repair target absent')
    a.output.mkdir(parents=True)
    descriptor = a.output / 'Log.txt'
    text = original.rstrip() + '\n\n'
    text += '# Constructed continuation descriptor; source run did not record these entries.\n'
    text += '# This records existing artifact locations, not completion of inference.\n'
    text += 'FN_Orthogroups: ' + str(clusters) + '\n'
    text += 'WorkingDirectory_Trees: ' + str(wd) + '/\n'
    descriptor.write_text(text)
    after = locate(a.output)
    if (after['clusters'] != str(clusters)
            or Path(after['tree_working_directory']) != wd
            or before['working_directory_chain'] != after['working_directory_chain']):
        raise ValueError('Native descriptor readback differs')
    unavailable = {key: path for key, path in after.items()
                   if key.endswith('species_tree') and not Path(path).is_file()}
    if sha(log) != log_hash:
        raise ValueError('Original log changed during preparation')
    result = dict(status='prepared_native_readback_continuation_descriptor',
                  source_results=str(source), original_log_sha256=log_hash,
                  universe_receipt_sha256=sha(a.universe_receipt),
                  cluster_sha256=sha(clusters), species_ids_sha256=sha(ids),
                  native_locator_before=before, native_locator_after=after,
                  missing_species_tree_artifacts=unavailable,
                  repair_target_empty_at_preparation=missing_tree.stat().st_size == 0,
                  parser_sha256=sha(Path(files.__file__)), script_sha256=sha(Path(__file__)),
                  artifacts={'Log.txt': sha(descriptor)},
                  next_requirements=['Completed and independently installed OG0000017 repair',
                                     'Reviewed rooted 526-taxon species tree supplied through -s, or verified species-tree artifacts',
                                     'Pinned native continuation command, resource plan and output protection',
                                     'Independent reconciliation and species-tree sensitivity validation'],
                  interpretation='Native read-only locator resolves audited clusters and existing gene-tree directory through separate explicitly constructed metadata. No completed species tree is asserted and no inference is launched. Original run log is unchanged.')
    (a.output / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
