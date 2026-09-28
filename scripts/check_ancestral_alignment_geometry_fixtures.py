#!/usr/bin/env python3
"""Compare full extant geometry distances for all frozen short-chain inputs."""
import argparse
from io import StringIO
import itertools
import json
from pathlib import Path
import subprocess
import tempfile
import time
from Bio import SeqIO
from ancestral_chain_attempt import sha, write_json
from ancestral_extant_alignment_geometry import project, signature, signature_distance
from readback_independent_baliphy_chain import blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    inputs = Path('results/ancestral/baliphy-independent-chain-inputs-20260927-v1/chain_inputs.json')
    proof = json.loads(Path('metadata/baliphy_sample_mapping_final_readback_completed_20260927.json').read_text())
    ap = Path(proof['audit_receipt'])
    assert sha(ap) == proof['audit_receipt_sha256']
    audit = json.loads(ap.read_text())
    selected = {}
    for row in json.loads(inputs.read_text()):
        selected.setdefault(row['effective_input_group'], row)
    assert len(selected) == 135
    binary = Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/bin/alignment-distances').resolve()
    pins = {str(p): sha(p) for p in [inputs, ap, binary, Path(__file__), Path('scripts/ancestral_extant_alignment_geometry.py')]}
    started = time.monotonic()
    results = []
    for group_id, row in sorted(selected.items()):
        root = Path('results/ancestral/baliphy-sample-mapping-20260927-v1') / row['original_configuration_ids'][0]
        rp = root / 'receipt.json'
        assert sha(rp) == audit['pins'][str(rp)]
        receipt = json.loads(rp.read_text())
        for name, digest in receipt['artifacts'].items():
            assert sha(root / name) == digest
        assert sha(row['alignment']) == row['alignment_sha256']
        observed = {r.id: str(r.seq).replace('-', '').upper() for r in SeqIO.parse(row['alignment'], 'fasta')}
        paths = list(root.glob('mapping-check-*/C1.P1.fastas'))
        assert len(paths) == 1
        encoded, iterations = [], []
        for iteration, text in blocks(paths[0]):
            records = list(SeqIO.parse(StringIO(text.lstrip()), 'fasta'))
            sequences = {r.id: str(r.seq).upper() for r in records}
            assert len(records) == len(sequences)
            encoded.append(project(sequences, observed)); iterations.append(iteration)
        assert iterations == [0, 10, 20]
        with tempfile.TemporaryDirectory() as temp:
            fasta = []
            for i, e in enumerate(encoded):
                p = Path(temp) / f'{i}.fasta'
                p.write_text(''.join(f'>{t}\n{s}\n' for t, s in e.items())); fasta.append(str(p))
            result = subprocess.run([str(binary), 'AxA', *fasta, '--alphabet', 'Amino-Acids', '--distances', 'pairwise'],
                                    capture_output=True, text=True, check=True, timeout=600)
        actual = [[float(v) for v in line.split()] for line in result.stdout.splitlines() if line.strip()]
        assert len(actual) == 3 and all(len(r) == 3 for r in actual)
        expected = [[0]*3 for _ in range(3)]
        for i in range(2):
            a = signature(encoded[i])
            for j in range(i+1,3):
                b = signature(encoded[j])
                expected[i][j] = expected[j][i] = signature_distance(a,b)
                del b
            del a
        rounded = [[float(format(v, '.6g')) for v in row] for row in expected]
        assert actual == rounded, (group_id, actual, expected, rounded)
        record = dict(group=group_id, input_receipt=str(rp), input_receipt_sha256=sha(rp),
            alignment_sha256=row['alignment_sha256'], sample_path=str(paths[0]), sample_sha256=sha(paths[0]),
            tips=len(observed), residues=sum(map(len,observed.values())), input_X=sum(s.count('X') for s in observed.values()),
            iterations=iterations, expected_exact_integer=expected, expected_six_significant_digits=rounded, actual_printed=actual, raw_stdout=result.stdout, stderr=result.stderr)
        results.append(record)
        write_json(args.output / f'group-{len(results):03d}.json', record)
        print(json.dumps(dict(groups=len(results), elapsed_seconds=time.monotonic()-started)), flush=True)
    for p, digest in pins.items():
        assert sha(p) == digest
    write_json(args.output / 'receipt.json', dict(status='passed_all_frozen_extant_geometry_distances_at_printed_precision',
        groups=len(results), alignments=3*len(results), matrix_entries=9*len(results), pins=pins,
        input_X_observations=sum(r['input_X']*3 for r in results), elapsed_seconds=time.monotonic()-started,
        artifacts={p.name:sha(p) for p in sorted(args.output.glob('group-*.json'))},
        scope='All 135 effective frozen short inputs, full extant tip sets, three saved alignments. '
        'A/gap encoding preserves positions including input X and discards amino-acid state information. '
        'Native comparison is at six significant printed digits, not exact integer equality. '
        'Not a production-chain convergence assessment or ancestral-state mixing diagnostic.'))


if __name__ == '__main__':
    main()
