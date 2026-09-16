#!/usr/bin/env python3
"""Independently verify full staged file coverage, bytes, and source isolation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from orthofinder.utils import util


def digest(p):
    with p.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    plan = json.loads(args.plan.read_text())
    base = Path(plan['output']).resolve()
    receipt = json.loads((base / 'receipt.json').read_text())
    if receipt['plan_sha256'] != digest(args.plan) or receipt['launch_ready']:
        raise ValueError('Unexpected staging receipt')
    with (Path(plan['inventory_receipt']).parent / 'families.tsv').open() as handle:
        trees = {r['family'] for r in csv.DictReader(handle, delimiter='\t') if r['tree_status'] == 'tip_branch_checks_passed'}
    source_wd = Path(plan['source_working_directories'][0])
    species = {int(line.split(':')[0]) for line in (source_wd / 'SpeciesIDs.txt').read_text().splitlines() if line.strip()}
    expected = {'WorkingDirectory/' + n for n in ['SpeciesIDs.txt','SequenceIDs.txt','clusters_OrthoFinder.txt_id_pairs.txt','SpeciesTree_unrooted_ids.txt']}
    expected |= {f'WorkingDirectory/Species{i}.fa' for i in species}
    expected |= {f'WorkingDirectory/Trees_ids/{og}.txt' for og in trees}
    expected = {'Source/' + n for n in expected} | {'species_tree_taxa.nwk'}
    reports = []
    for report in receipt['guides']:
        folder = base / report['guide']
        manifest = folder / 'copied_files.tsv'
        if digest(manifest) != report['manifest_sha256']:
            raise ValueError('Changed manifest')
        with manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        names = {r['relative_path'] for r in rows}
        if names != expected or len(names) != len(rows):
            raise ValueError('Incomplete or duplicated file coverage')
        actual = {str(p.relative_to(folder)) for p in (folder / 'Source').rglob('*') if p.is_file()}
        if actual != expected - {'species_tree_taxa.nwk'} | {'Source/Log.pending.txt'}:
            raise ValueError('Unexpected files in staged source')
        total = 0
        for row in rows:
            copy = folder / row['relative_path']
            original = Path(row['source'])
            if copy.is_symlink() or (copy.stat().st_dev,copy.stat().st_ino) == (original.stat().st_dev,original.stat().st_ino):
                raise ValueError('Source alias')
            if digest(copy) != row['sha256'] or digest(original) != row['sha256'] or copy.stat().st_size != int(row['bytes']):
                raise ValueError('Source/copy digest or size mismatch')
            total += copy.stat().st_size
        if total != report['copied_bytes'] or len(rows) != report['copied_files']:
            raise ValueError('Summary mismatch')
        # Native sequence accounting verifies every staged species FASTA is found.
        info = util.GetSeqsInfo([str(folder / 'Source/WorkingDirectory')], sorted(species), max(species)+1)
        if set(info.nSeqsPerSpecies) != species or any(n <= 0 for n in info.nSeqsPerSpecies.values()):
            raise ValueError('Native species sequence coverage differs')
        reports.append(dict(guide=report['guide'],checked_files=len(rows),checked_bytes=total,source_aliases=0,
                            native_species_fasta_count=len(info.nSeqsPerSpecies),native_fasta_sequences=info.nSeqs))
    result = dict(status='passed_complete_staged_file_readback',guides=reports,
                  receipt_sha256=digest(base / 'receipt.json'),script_sha256=digest(Path(__file__)),
                  launch_ready=False,pending_families=receipt['pending_families'],
                  scope='All copied files independently hashed against original files; exact expected file universe and independent inode identity checked. Validates staging, not orthology or native runtime access completeness.')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
