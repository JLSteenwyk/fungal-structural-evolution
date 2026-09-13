#!/usr/bin/env python3
"""Fit sequence and native-3Di branches on matched marker sequence topologies."""
import argparse
import csv
import fcntl
import hashlib
import json
import math
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import Phylo, SeqIO
from compare_marker_structures import sha
from assess_pae_sensitivity import checked_receipt


def tree_edges(path, taxa):
    tree = Phylo.read(path, 'newick')
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != len(taxa) or set(tips) != taxa:
        raise ValueError('Tree tips differ from input')
    result = {}
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        side = {t.name for t in clade.get_terminals()}
        a, b = tuple(sorted(side)), tuple(sorted(taxa - side))
        key = min((a, b), key=lambda x: (len(x), x))
        if clade.branch_length is None or not math.isfinite(clade.branch_length) or clade.branch_length < 0:
            raise ValueError('Invalid branch length')
        # Suppress a degree-two display root by summing its two halves.
        result[key] = result.get(key, 0.) + clade.branch_length
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'models', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.inputs, args.models, args.output = args.inputs.resolve(), args.models.resolve(), args.output.resolve()
    checked_receipt(args.inputs)
    checked_receipt(args.models)
    executable = shutil.which('iqtree3')
    if executable is None:
        raise FileNotFoundError('iqtree3')
    version = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True).stdout
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / '.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = {'input_receipt_sha256': sha(args.inputs / 'receipt.json'),
        'model_receipt_sha256': sha(args.models / 'receipt.json'),
        'script_sha256': sha(Path(__file__)), 'executable': executable,
        'executable_sha256': sha(Path(executable)), 'version': version,
        'workers': 4, 'threads_per_fit': 1, 'memory_per_fit': '2G',
        'planning': '52 available markers, 4-12 taxa, 91-1181 columns, four fits each; provisional 0.1-10 minutes per fit, 0.1-9 hours wall time at four workers, 2 GB output headroom. Local CPU only, no paid services.',
        'analysis': 'AA LG+F+G4 tree search; AF+G4, AF+F+G4 and LLM+G4 structural fits on that fixed unrooted AA topology. Identical sequences retained. Branches are expected state substitutions/site, not Angstrom or change/year. No branch uncertainty yet.'}
    config_path = args.output / 'config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError('Changed run configuration; use a new output')
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    marker_rows = list(csv.DictReader((args.inputs / 'marker_summary.tsv').open(), delimiter='\t'))

    def run(row):
        marker = row['marker']
        folder = args.output / marker
        folder.mkdir(exist_ok=True)
        aa_path = args.inputs / marker / 'aa.faa'
        state_path = args.inputs / marker / '3di.faa'
        taxa = {r.id for r in SeqIO.parse(aa_path, 'fasta')}
        definitions = [('aa', aa_path, 'LG+F+G4'),
                       ('3di_af', state_path, str(args.models / 'Q.3Di.AF') + '+G4'),
                       ('3di_af_empirical', state_path, str(args.models / 'Q.3Di.AF') + '+F+G4'),
                       ('3di_llm', state_path, str(args.models / 'Q.3Di.LLM') + '+G4')]
        outputs = {}
        reference = None
        for label, alignment, model in definitions:
            prefix = folder / label
            tree_path = folder / (label + '.treefile')
            seed = int(hashlib.sha256(marker.encode()).hexdigest()[:8], 16) % 2147483646 + 1
            command = [executable, '-s', str(alignment), '-st', 'AA', '-m', model,
                       '-T', '1', '--mem', '2G', '--seed', str(seed), '-keep-ident', '--prefix', str(prefix)]
            if label != 'aa':
                command += ['-te', str(folder / 'aa.treefile')]
            run_config = {'command': command, 'alignment_sha256': sha(alignment),
                          'parent_config_sha256': sha(config_path),
                          'topology_sha256': sha(folder / 'aa.treefile') if label != 'aa' else None}
            run_path = folder / (label + '.config.json')
            if run_path.exists() and json.loads(run_path.read_text()) != run_config:
                raise ValueError('Changed fit inputs or command')
            run_path.write_text(json.dumps(run_config, indent=2) + '\n')
            receipt_path = folder / (label + '.receipt.json')
            if receipt_path.exists():
                result = json.loads(receipt_path.read_text())
                if result['config_sha256'] != sha(run_path):
                    raise ValueError('Changed completed fit configuration')
                for name, checksum in result['artifacts'].items():
                    if sha(folder / name) != checksum:
                        raise ValueError('Changed completed fit output')
            else:
                start = time.monotonic()
                with (folder / (label + '.stdout.log')).open('a') as log:
                    process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
                if process.returncode != 0:
                    raise RuntimeError(f'{marker}/{label} failed with exit {process.returncode}; inspect log')
                tree_edges(tree_path, taxa)
                result = {'status': 'completed_point_estimate', 'elapsed_seconds': time.monotonic() - start,
                    'config_sha256': sha(run_path), 'artifacts': {label + suffix: sha(folder / (label + suffix))
                    for suffix in ['.treefile', '.iqtree', '.log']}}
                receipt_path.write_text(json.dumps(result, indent=2) + '\n')
            edges = tree_edges(tree_path, taxa)
            if reference is None:
                reference = set(edges)
            elif set(edges) != reference:
                raise ValueError('Fixed sequence topology changed during structural fit')
            outputs[label] = edges
        records = []
        taxon_set_hash = hashlib.sha256('\n'.join(sorted(taxa)).encode()).hexdigest()
        for split in sorted(reference):
            records.append({'marker': marker, 'taxon_set_sha256': taxon_set_hash,
                'split_taxa': ','.join(split), 'split_taxon_count': len(split),
                **{label + '_branch_length': edges[split] for label, edges in outputs.items()}})
        with (folder / 'paired_branches.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(records)
        print(marker, 'completed', len(taxa), len(records), flush=True)
        return {'marker': marker, 'taxa': len(taxa), 'branches': len(records),
            'paired_branches_sha256': sha(folder / 'paired_branches.tsv')}

    results = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run, row) for row in marker_rows if row['status'] == 'ready_for_inference']
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except Exception:
            for future in futures:
                future.cancel()
            raise
    receipt = {'status': 'complete_matched_topology_point_estimates', 'config_sha256': sha(config_path),
        'markers': len(results), 'results': sorted(results, key=lambda r: r['marker']),
        'interpretation': 'Point estimates conditional on sequence-derived marker topology and confidence-selected sites. Shared ancestry, sparse coverage, correlated 3Di features, support/model adequacy and branch uncertainty prevent acceleration claims at this stage.'}
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('Completed all', len(results), 'markers', flush=True)


if __name__ == '__main__':
    main()
