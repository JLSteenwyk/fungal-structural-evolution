#!/usr/bin/env python3
"""Recreate the tied-quantile fixture from the preserved raw-log failure record."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from run_ortholog_pair_guide_comparison import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--failure', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); record = json.loads(args.failure.read_text())['diagnosis']
    for path, digest in record['source_pins'].items():
        assert sha(path) == digest
    output = record['report']; assert sha(output['path']) == output['sha256']
    report = json.loads(Path(output['path']).read_text())
    assert report['pins'] == record['source_pins']
    assert report['variables'][record['variable']] == record['original']
    paths = [p for p in report['pins'] if Path(p).name == 'manifest-' + record['cutoff'] + '.json']
    assert len(paths) == 1
    manifest = json.loads(Path(paths[0]).read_text()); values = []
    assert len(manifest['chains']) == len({r['seed'] for r in manifest['chains']}) == 4
    for chain in manifest['chains']:
        assert chain['model_input_identity'] == record['group']
        with Path(chain['log']).open() as f:
            reader = csv.reader(f, delimiter='\t'); header = next(reader)
            rows = list(reader)
        iteration = header.index('iter'); variable = header.index(record['variable'])
        assert [int(r[iteration]) for r in rows] == manifest['expected_iterations']
        values.append([float(r[variable]) for r in rows if int(r[iteration]) > int(record['cutoff'])])
    assert np.asarray(values).shape == (4, 750)
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.output / 'trace.npy'; np.save(path, np.asarray(values, dtype=float))
    assert sha(path) == record['trace_sha256']
    source = dict(record, trace_path=str(path))
    (args.output / 'source.json').write_text(json.dumps(source, indent=2) + '\n')
    print(json.dumps(dict(status='recreated_exact_source_tied_quantile_regression', path=str(path), sha256=sha(path))))


if __name__ == '__main__':
    main()
