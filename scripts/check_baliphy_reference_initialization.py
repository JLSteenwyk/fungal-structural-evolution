#!/usr/bin/env python3
"""Native software fixtures for reference startup; no fungal pilot or inference."""
import argparse
import copy
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path
import re
import subprocess
import time

from Bio import SeqIO

from ancestral_chain_attempt import run_attempt, sha, write_json
from baliphy_reference_initialization import alignment_records, homology_signature, transform, validate_initial_json


def rejects(action):
    try:
        action()
    except (AssertionError, KeyError):
        return
    raise AssertionError('Invalid fixture was accepted')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    root = a.output.resolve(); root.mkdir(exist_ok=False)
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    cg = Path('/sys/fs/cgroup') / group.lstrip('/')
    limits = {k: (cg / k).read_text().strip() for k in ['cpu.max', 'memory.max', 'memory.swap.max']}
    assert limits == {'cpu.max': '200000 100000', 'memory.max': str(16 * 2**30), 'memory.swap.max': '0'}
    assert all(os.environ.get(k) == '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'])
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/bali-phy').resolve()
    chain_path = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    chains = json.loads(chain_path.read_text())
    alignment = root / 'alignment.faa'; tree = root / 'tree.nwk'
    alignment.write_text('>a\nA-CDE-FG\n>b\nA-C-EYFG\n>c\nATCDE-F-\n>d\n--CDEF-G\n')
    tree.write_text('((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);\n')
    reference = alignment_records(alignment); native = []; negatives = []; bindings = {}
    started = time.monotonic()
    def execute(program, seed, name, test=True):
        command = ['/usr/bin/prlimit', '--as=' + str(8 * 2**30), '--cpu=300',
                   '--fsize=' + str(64 * 2**20), '--', str(binary), '--seed', str(seed),
                   'run', str(program), '--log-format', 'json']
        command += ['--test'] if test else ['--iterations', '20', '--name', 'kernel-check']
        config = dict(command=command, timeout_seconds=360,
                      pins={str(x): sha(x) for x in [Path('/usr/bin/prlimit'), binary, program, alignment, tree]})
        receipt = run_attempt(root / name, config)
        r = json.loads(receipt.read_text()); assert r['exit_code'] == 0, receipt
        bindings[str(receipt)] = sha(receipt)
        for x, h in r['artifacts'].items(): bindings[str(receipt.parent / x)] = h
        if test:
            return json.loads((receipt.parent / 'stdout.log').read_text())
        return receipt.parent
    for prior in ['broad', 'centered', 'package']:
        c = next(x for x in chains if x['prior_label'] == prior)
        original = Path(c['program']); assert sha(original) == c['program_sha256']
        bindings[str(original)] = sha(original)
        old = original.read_text().replace(str(Path(c['alignment']).resolve()), str(alignment)).replace(str(Path(c['tree']).resolve()), str(tree))
        original_fixture = root / (prior + '-original.hs'); original_fixture.write_text(old)
        modified_fixture = root / (prior + '-reference.hs'); modified_fixture.write_text(transform(old))
        for role in range(1, 5):
            seed = 425 + role
            baseline = execute(original_fixture, seed, prior + '-baseline-' + str(role))
            value = execute(modified_fixture, seed, prior + '-reference-' + str(role))
            result = validate_initial_json(value, reference, ['a', 'b', 'c', 'd'])
            for key in ['S1/', 'I1/']:
                assert value['parameters/'][key] == baseline['parameters/'][key], 'Initial parameter draw changed'
            native.append(dict(prior_label=prior, role=role, seed=seed, same_seed_parameters_preserved=True, **result))
            print('native_software_initialization', prior, role, flush=True)
        folder = execute(modified_fixture, 426, prior + '-kernel', test=False)
        fasta = list(folder.glob('kernel-check-*/C1.P1.fastas')); assert len(fasta) == 1
        blocks = re.split(r'^iterations = (\d+)\s*\n', fasta[0].read_text(), flags=re.M)
        assert [int(blocks[i]) for i in range(1, len(blocks), 2)] == [0, 10, 20]
        geometries = []
        for i in range(1, len(blocks), 2):
            seqs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(io.StringIO(blocks[i + 1]), 'fasta')]
            extant = [(label, s) for label, s in seqs if label in ['a', 'b', 'c', 'd']]
            assert len(extant) == 4
            assert {label: s.replace('-', '') for label, s in extant} == {label: s.replace('-', '') for label, s in reference}
            geometry = homology_signature(extant); geometries.append(geometry)
            if i == 1: assert geometry == homology_signature(reference)
        assert any(g != geometries[0] for g in geometries[1:]), 'Alignment transition did not move in software fixture'
        logged = [json.loads(x) for x in (fasta[0].parent / 'C1.log.json').read_text().splitlines()]
        assert logged[0] == dict(fields=['iter', 'prior', 'likelihood', 'posterior'],
                                 nested=True, format='MCON', version='0.2')
        logs = logged[1:]
        assert len(logs) == 21 and [x['iter'] for x in logs] == list(range(21))
        assert all('referenceInitialization/' not in x['parameters//'] for x in logs), 'Test-only audit leaked into inference logs'
        native[-1]['kernel_fixture_twenty_iterations_changed_alignment'] = True
        for name, change in [
            ('ancestral_lengths_fixed', lambda x: x['parameters/']['referenceInitialization/'].update(fixedLengthNodes=x['parameters/']['referenceInitialization/']['allLengthNodes'])),
            ('density_mismatch', lambda x: x['parameters/']['referenceInitialization/'].update(priorFromDirectFormula=0)),
            ('score_mismatch', lambda x: x.update(posterior=x['posterior'] + 10)),
            ('reference_geometry_mismatch', lambda x: x['parameters/']['referenceInitialization/'].update(extantAlignment='>a\nACDEFG--\n>b\nACEYFG--\n>c\nATCDEF--\n>d\nCDEFG---\n')),
        ]:
            bad = copy.deepcopy(value); change(bad)
            rejects(lambda: validate_initial_json(bad, reference, ['a', 'b', 'c', 'd']))
            negatives.append(prior + ':' + name)
    for prior in ['broad', 'centered', 'package']:
        source = (root / (prior + '-original.hs')).read_text()
        for before in ['import Probability.Random\n', ';return (parameterLogValues loggerValues)', ';(alignment,properties_A) <- sampleWithProps (phyloAlignment tree imodel scale1 sequence_lengths)']:
            rejects(lambda: transform(source.replace(before, '')))
            negatives.append(prior + ':missing_transform_anchor:' + before[:24])
    for name, content in [('duplicate_labels', '>a\nAC\n>a\nAC\n'), ('ragged_matrix', '>a\nAC\n>b\nA\n')]:
        path = root / (name + '.faa'); path.write_text(content)
        rejects(lambda: alignment_records(path)); negatives.append(name)
    for x, h in bindings.items(): assert sha(x) == h, x
    paths = [Path(__file__), Path('scripts/baliphy_reference_initialization.py'),
             Path('scripts/ancestral_chain_attempt.py'), chain_path, binary, Path('/usr/bin/prlimit')]
    bindings.update({str(x): sha(x) for x in paths})
    receipt = dict(status='passed_native_reference_initialization_software_contracts',
                   checked_utc=datetime.now(timezone.utc).isoformat(), native_initializations=12,
                   original_initializations=12, kernel_fixtures=3, iterations_per_kernel_fixture=20,
                   native_results=native, invalid_inputs_and_audits_rejected=negatives,
                   source_hashes=bindings, elapsed_seconds=time.monotonic() - started,
                   cgroup_limits=limits, binary_sha256=sha(binary),
                   scientific_eligibility=False,
                   scope='Synthetic four-tip software fixtures only. Same-seed stochastic parameter draws, extant homology, fixed-tip/free-ancestor lengths, native degree-aware alignment density, and moving alignment kernels checked under caps. No fungal pilot, new posterior, convergence or memory repair claim.')
    with a.receipt.open('x') as out: out.write(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ['source_hashes', 'native_results']}), flush=True)


if __name__ == '__main__': main()
