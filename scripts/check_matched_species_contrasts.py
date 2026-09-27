#!/usr/bin/env python3
"""Check every species contrast from raw record keys and exported endpoint taxa."""
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy import sparse
from screen_duplication_domain_alignment_coverage import sha


def main():
    root = Path('results/phylogeny/matched-species-contrasts-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    for path, digest in receipt['pins'].items():
        assert sha(path) == digest, path
    for name, digest in receipt['artifacts'].items():
        assert sha(root / name) == digest, name
    taxa = json.loads((root / 'taxa.json').read_text())
    index = {name: i for i,name in enumerate(taxa)}
    nodes = {}
    with open('results/orthology/selected-taxon-inputs-20260927-v1/selected_node_taxa.tsv') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            nodes[row['role'], row['node_id']] = row
    counts = Counter()
    with gzip.open('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            counts[row['target_id'], row['background_id']] += 1
    patterns = {}
    matrix = sparse.load_npz(root / 'species_contrast_design.npz')
    for row in csv.DictReader((root / 'patterns.tsv').open(), delimiter='\t'):
        indices = [index[t] for t in json.loads(row['taxa'])]
        weights = json.loads(row['twice_weights'])
        assert indices == sorted(set(indices)) and all(w != 0 for w in weights)
        raw = list(zip(indices, weights))
        label = hashlib.sha256(json.dumps(raw, separators=(',', ':')).encode()).hexdigest()
        assert label == row['species_pattern_id'] and label not in patterns
        r = int(row['row_index'])
        vector = matrix.getrow(r)
        np.testing.assert_array_equal(vector.indices, indices)
        np.testing.assert_array_equal(vector.data * 2, weights)
        assert sum(weights) == 0
        patterns[label] = (r, dict(zip(indices, weights)))
    assert sorted(r for r,w in patterns.values()) == list(range(matrix.shape[0]))
    seen = set()
    used_patterns = set()
    with gzip.open(root / 'selected_pair_patterns.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = row['target_id'], row['background_id']
            assert key not in seen; seen.add(key)
            assert int(row['selected_record_uses']) == counts[key]
            t, b = nodes['target',key[0]], nodes['background',key[1]]
            assert t['taxon_a'] == t['taxon_b']
            assert (t['guide'], t['family']) == (b['guide'], b['family'])
            a, left, right = index[t['taxon_a']], index[b['taxon_a']], index[b['taxon_b']]
            assert [a,left,right] == [int(row[k]) for k in ['focal_index','background_index_a','background_index_b']]
            expected = Counter({a: 2})
            expected[left] -= 1; expected[right] -= 1
            expected = {i:w for i,w in expected.items() if w}
            label = row['species_pattern_id']; used_patterns.add(label)
            assert patterns[label][1] == expected
    assert seen == set(counts) and used_patterns == set(patterns)
    assert len(seen) == receipt['unique_selected_pairs'] and sum(counts.values()) == receipt['selected_records']
    kernels = Path('results/phylogeny/species-distance-kernels-20260927-v1')
    distances = {path.stem: np.load(path)['distances'] for path in kernels.glob('*.npz')}
    observed = set()
    for row in csv.DictReader((root / 'kernel_quadratic_forms.tsv').open(), delimiter='\t'):
        key = row['tree'], row['species_pattern_id']
        assert key not in observed; observed.add(key)
        weights = patterns[key[1]][1]
        # General zero-sum expression, without the producer's special three-tip formula.
        expected = -.125 * sum(w * v * distances[key[0]][i,j] for i,w in weights.items() for j,v in weights.items())
        np.testing.assert_allclose(float(row['kernel_quadratic_form']), expected, rtol=1e-10, atol=1e-12)
    assert observed == {(tree,label) for tree in distances for label in patterns}
    active = np.flatnonzero(np.asarray(abs(matrix).sum(axis=0)).ravel())
    singular = np.linalg.svd(matrix[:, active].toarray(), compute_uv=False)
    tolerance = max(matrix.shape) * np.finfo(float).eps * singular[0]
    rank = int((singular > tolerance).sum())
    assert rank <= len(active) - 1
    result = dict(status='passed_full_matched_species_contrast_readback',
                  selected_records=sum(counts.values()), selected_pairs=len(seen), patterns=len(patterns),
                  quadratic_forms_checked=len(observed), active_taxa=len(active), contrast_rank=rank,
                  svd_rank_tolerance=float(tolerance), singular_values=singular.tolist(),
                  source_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(__file__),
                  scope='Every original record multiplicity, taxon weight, sparse row and quadratic form checked. Rank characterizes the full pre-qualification working design, not sample size, model adequacy or a fitted evolutionary effect.')
    Path('metadata/matched_species_contrast_readback_20260927.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'singular_values'}, indent=2))


if __name__ == '__main__':
    main()
