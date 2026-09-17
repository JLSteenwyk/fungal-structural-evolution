#!/usr/bin/env python3
"""Exercise the production family runner on synthetic sequences and invalid outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import random
from run_expanded_family_trees import run_task, check_alignment, check_tree, sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    binary = Path('.cache/envs/orthofinder/lib/python3.12/site-packages/orthofinder/bin').resolve()
    plan = dict(output=str(a.output / 'jobs'), inputs=str(a.output), famsa=str(binary / 'famsa'),
                fasttree=str(binary / 'FastTree'), resources=dict(step_timeout_hours=1))
    rng = random.Random(20260917)
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    results = []
    for n in [3, 12]:
        ancestor = ''.join(rng.choice(alphabet) for _ in range(600))
        source = a.output / f'synthetic{n}.faa'
        with source.open('w') as handle:
            for i in range(n):
                seq = list(ancestor)
                # Two identical sequences test tip preservation through deduplication.
                if i >= 2:
                    for j in rng.sample(range(600), 20):
                        seq[j] = rng.choice(alphabet)
                if n == 12 and i == 11:
                    seq += [rng.choice(alphabet) for _ in range(300)]
                handle.write(f'>0_{i}\n' + ''.join(seq) + '\n')
        key = hashlib.sha256(('\n'.join(sorted(f'0_{i}' for i in range(n))) + '\n').encode()).hexdigest()
        row = dict(membership_sha256=key, relative_path=source.name, sha256=sha(source), proteins=n)
        result = run_task(row, plan, 'synthetic-software-fixture-v1')
        assert run_task(row, plan, 'synthetic-software-fixture-v1') == result
        if n == 12:
            assert result['dimensions']['trimmed_columns'] < result['dimensions']['raw_columns']
        folder = a.output / 'jobs' / 'families' / key
        badtree = a.output / f'bad{n}.nwk'
        badtree.write_text('(wrong:0.1,0_1:0.2,0_2:0.3);\n')
        try:
            check_tree(badtree, source)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid tree accepted')
        badraw = a.output / f'bad{n}.faa'
        rawtext = (folder / 'raw.faa').read_text()
        badraw.write_text(rawtext.replace('>0_0', '>wrong', 1))
        try:
            check_alignment(source, badraw, folder / 'trimmed.faa')
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid alignment accepted')
        results.append(dict(proteins=n, dimensions=result['dimensions'], checkpoint_reuse_passed=True,
                            wrong_tree_tips_rejected=True, wrong_alignment_ids_rejected=True))
    receipt = dict(status='passed_synthetic_family_runner_checks', cases=results,
                   runner_sha256=sha(Path('scripts/run_expanded_family_trees.py')),
                   script_sha256=sha(Path(__file__)),
                   scope='Software validation only, including duplicate tips, actual trimming and checkpoint reuse; not a biological pilot or full-scale runtime benchmark.')
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
