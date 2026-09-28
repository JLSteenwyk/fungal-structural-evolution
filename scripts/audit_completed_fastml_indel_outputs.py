#!/usr/bin/env python3
"""Audit terminal FastML outputs with explicit root-comment handling.

Do not change failed producer receipts or claim independent likelihood replay.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
from Bio import Phylo, SeqIO
from ancestral_chain_attempt import sha, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    assert not args.output.exists()
    unit = 'fungal-ancestral-fastml-indels-20260927.service'
    state = dict(line.split('=', 1) for line in subprocess.check_output(
        ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'MainPID',
         '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
    assert state == dict(ActiveState='failed', MainPID='0', Result='exit-code', ExecMainStatus='1')
    pp = Path('metadata/ancestral_fastml_indel_plan_20260927.json')
    plan = json.loads(pp.read_text())
    for path, digest in plan['pins'].items():
        assert sha(path) == digest
    root = Path(plan['output'])
    terminal = json.loads((root / 'receipt.json').read_text())
    assert terminal['failed_jobs'] == 153 and len(terminal['jobs']) == 156
    jobs = []
    for job in plan['jobs']:
        folder = root / job['job_id']
        rp = folder / 'receipt.json'
        assert sha(rp) == terminal['job_receipts'][str(rp.relative_to(root))]
        receipt = json.loads(rp.read_text())
        assert receipt['job'] == job and receipt['plan_sha256'] == sha(pp)
        for name, digest in receipt['artifacts'].items():
            assert sha(folder / name) == digest
        row = dict(job_id=job['job_id'], producer_status=receipt['status'],
                   source_receipt_sha256=sha(rp))
        if not job['character_count']:
            assert receipt['status'] == 'no_coded_characters_no_indel_inference'
            row['readback_status'] = 'no_coded_characters'
            jobs.append(row)
            continue
        assert receipt['exit_code'] == 0
        tree = Phylo.read(folder / 'RESULTS/TheTree.INodes.ph', 'newick')
        assert tree.root.name is None and tree.root.comment == 'N1'
        tree.root.name = tree.root.comment
        nodes = [n.name for n in tree.find_clades()]
        assert None not in nodes and len(nodes) == len(set(nodes))
        seq = {r.id: str(r.seq) for r in SeqIO.parse(job['characters'], 'fasta')}
        assert len(seq) == job['proteins']
        assert {n.name for n in tree.get_terminals()} == set(seq)
        assert {len(s) for s in seq.values()} == {job['character_count']}
        assert all(set(s) <= set('01?') for s in seq.values())
        seen, maximum, known = set(), 0., 0
        with (folder / 'RESULTS/AncestralReconstructPosterior.txt').open() as handle:
            for r in csv.DictReader(handle, delimiter='\t'):
                pos, node, value = int(r['POS']), r['Node'], float(r['Prob'])
                assert r['State'] == '1' and math.isfinite(value) and 0 <= value <= 1
                assert node in nodes and 1 <= pos <= job['character_count']
                assert (pos, node) not in seen
                seen.add((pos, node))
                if node in seq and seq[node][pos - 1] != '?':
                    maximum = max(maximum, abs(value - int(seq[node][pos - 1])))
                    known += 1
        assert len(seen) == len(nodes) * job['character_count']
        assert maximum < 1e-6
        row.update(readback_status='passed_output_identity_bounds_grid_and_known_tip_checks',
                   root_label_source='Newick root comment N1', probability_rows=len(seen),
                   known_tip_probabilities=known, maximum_known_tip_error=maximum)
        jobs.append(row)
    assert len(jobs) == 156
    write_json(args.output, dict(status='completed_output_readback_not_likelihood_validation',
        producer_terminal_state=state, plan_sha256=sha(pp),
        terminal_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(__file__),
        jobs=jobs, probability_rows=sum(r.get('probability_rows', 0) for r in jobs),
        known_tip_probabilities=sum(r.get('known_tip_probabilities', 0) for r in jobs),
        scope='Original 153 validator failures retained. Root comment N1 explicitly restored only '
              'in readback; identity, probability bounds, complete grids and known tips checked. '
              'No independent likelihood or internal probability replay, optimization qualification, '
              'ascertainment correction resolution, or posterior ensemble qualification.'))
    print('Checked', len(jobs), 'inputs;', sum(r.get('probability_rows', 0) for r in jobs), 'probability rows')


if __name__ == '__main__':
    main()
