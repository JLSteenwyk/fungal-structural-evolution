#!/usr/bin/env python3
"""Read finished FastML output, including its root label stored as a comment."""
import csv
import json
import math
from pathlib import Path
from Bio import Phylo, SeqIO
from prepare_case_ancestral_neighborhoods import sha


def main():
    plan_path = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    plan = json.loads(plan_path.read_text())
    root = Path(plan['output'])
    rows = []
    for rp in sorted(root.glob('*/receipt.json')):
        receipt = json.loads(rp.read_text())
        assert receipt['plan_sha256'] == sha(plan_path)
        job = receipt['job']
        for relative, digest in receipt['artifacts'].items():
            assert sha(rp.parent / relative) == digest
        if not job['character_count']:
            rows.append(dict(job_id=job['job_id'], status='no_coded_characters'))
            continue
        assert receipt['exit_code'] == 0
        tree = Phylo.read(rp.parent / 'RESULTS/TheTree.INodes.ph', 'newick')
        assert tree.root.name is None and tree.root.comment == 'N1'
        tree.root.name = tree.root.comment
        names = [n.name for n in tree.find_clades()]
        assert None not in names and len(names) == len(set(names))
        seq = {r.id: str(r.seq) for r in SeqIO.parse(job['characters'], 'fasta')}
        assert {n.name for n in tree.get_terminals()} == set(seq)
        seen, tip_error = set(), 0.
        with (rp.parent / 'RESULTS/AncestralReconstructPosterior.txt').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                pos, node, probability = int(row['POS']), row['Node'], float(row['Prob'])
                assert row['State'] == '1' and math.isfinite(probability) and 0 <= probability <= 1
                assert 1 <= pos <= job['character_count'] and node in names
                assert (pos, node) not in seen
                seen.add((pos, node))
                if node in seq and seq[node][pos-1] != '?':
                    tip_error = max(tip_error, abs(probability-int(seq[node][pos-1])))
        assert len(seen) == len(names)*job['character_count'] and tip_error < 1e-6
        rows.append(dict(job_id=job['job_id'], status='output_shape_and_tip_states_verified',
                         probability_rows=len(seen), maximum_observed_tip_error=tip_error,
                         source_receipt_sha256=sha(rp)))
    output = Path('metadata/ancestral_fastml_indel_output_readback_20260927.json')
    output.write_text(json.dumps(dict(status='output_readback_snapshot', jobs=rows,
        inputs_with_receipts=len(rows), planned_inputs=len(plan['jobs']),
        plan_sha256=sha(plan_path), script_sha256=sha(__file__),
        scope='Root comment normalized for reading only. Original outputs and failed producer checks preserved. No independent internal-node probability/likelihood verification or terminal process success claim.'), indent=2)+'\n')
    print(len(rows), 'finished inputs read back')


if __name__ == '__main__':
    main()
