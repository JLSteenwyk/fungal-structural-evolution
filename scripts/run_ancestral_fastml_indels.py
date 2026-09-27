#!/usr/bin/env python3
"""Run all SIC matrices under the named FastML all-taxa correction sensitivity."""
import csv
import json
import math
import subprocess
from pathlib import Path
from Bio import Phylo, SeqIO
from prepare_case_ancestral_neighborhoods import sha


def main():
    pp = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    plan = json.loads(pp.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    out = Path(plan['output']).resolve()
    out.mkdir(parents=True, exist_ok=False)
    summaries = []
    for job in plan['jobs']:
        folder = out / job['job_id']
        folder.mkdir()
        seq = {r.id: str(r.seq) for r in SeqIO.parse(job['characters'], 'fasta')}
        assert len(seq) == job['proteins']
        assert {len(s) for s in seq.values()} == {job['character_count']}
        tree = Phylo.read(job['tree'], 'newick')
        assert {n.name for n in tree.get_terminals()} == set(seq)
        for node in tree.get_nonterminals():
            node.name = None
            node.confidence = None
        tree_path = folder / 'input_tree.nwk'
        Phylo.write(tree, tree_path, 'newick', format_branch_length='%1.12g')
        result = dict(job=job, plan_sha256=sha(pp))
        if not job['character_count']:
            result['status'] = 'no_coded_characters_no_indel_inference'
        else:
            settings = dict(plan['settings'])
            settings.update(_seqFile=str(Path(job['characters']).resolve()),
                            _treeFile=str(tree_path), _outDir=str(folder / 'RESULTS'))
            params = folder / 'parameters.txt'
            params.write_text(''.join(f'{k} {v}\n' for k,v in settings.items()))
            with (folder / 'stdout.log').open('w') as handle:
                run = subprocess.run([plan['binary'], str(params)], cwd=folder,
                                     stdout=handle, stderr=subprocess.STDOUT)
            result['exit_code'] = run.returncode
            try:
                assert run.returncode == 0, f'gainLoss exit {run.returncode}'
                fitted = Phylo.read(folder / 'RESULTS/TheTree.INodes.ph', 'newick')
                assert {n.name for n in fitted.get_terminals()} == set(seq)
                nodes = [n.name for n in fitted.find_clades()]
                assert None not in nodes and len(set(nodes)) == len(nodes)
                seen = set()
                max_tip_error = 0.
                with (folder / 'RESULTS/AncestralReconstructPosterior.txt').open() as handle:
                    for row in csv.DictReader(handle, delimiter='\t'):
                        pos, node, prob = int(row['POS']), row['Node'], float(row['Prob'])
                        assert row['State'] == '1' and math.isfinite(prob) and 0 <= prob <= 1
                        assert 1 <= pos <= job['character_count'] and node in nodes
                        key = pos, node
                        assert key not in seen
                        seen.add(key)
                        if node in seq and seq[node][pos-1] != '?':
                            error = abs(prob - int(seq[node][pos-1]))
                            max_tip_error = max(max_tip_error, error)
                assert len(seen) == job['character_count'] * len(nodes)
                assert max_tip_error < 1e-6, max_tip_error
                result.update(status='produced_all_node_probabilities_pending_independent_audit',
                              probability_rows=len(seen), maximum_observed_tip_error=max_tip_error)
            except Exception as error:
                result.update(status='failed_output_validation', error=repr(error))
        result['artifacts'] = {str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file()}
        (folder / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        summaries.append(dict(job_id=job['job_id'], status=result['status']))
        print(job['job_id'], result['status'], result.get('error', ''), flush=True)
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    receipt = dict(status='all_156_inputs_attempted_pending_independent_scientific_audit',
                   jobs=summaries, plan_sha256=sha(pp),
                   failed_jobs=sum(r['status']=='failed_output_validation' for r in summaries),
                   job_receipts={str(p.relative_to(out)):sha(p) for p in out.glob('*/receipt.json')},
                   scope=plan['scope'])
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    assert receipt['failed_jobs'] == 0, receipt['failed_jobs']


if __name__ == '__main__':
    main()
