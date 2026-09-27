#!/usr/bin/env python3
"""Export explicit additive endpoint species contrasts, without fitting a model."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from screen_duplication_domain_alignment_coverage import sha


def main():
    source = Path('results/orthology/selected-taxon-inputs-20260927-v1')
    kernels = Path('results/phylogeny/species-distance-kernels-20260927-v1')
    selections = Path('results/orthology/background-control-selection-20260927-v1/selections.tsv.gz')
    out = Path('results/phylogeny/matched-species-contrasts-20260927-v1')
    pins = {}
    for root in [source, kernels]:
        pins[str(root / 'receipt.json')] = sha(root / 'receipt.json')
        for name, digest in json.loads((root / 'receipt.json').read_text())['artifacts'].items():
            assert sha(root / name) == digest
            pins[str(root / name)] = digest
    source_receipt = json.loads((source / 'receipt.json').read_text())
    assert sha(selections) == source_receipt['pins'][str(selections)]
    pins[str(selections)] = sha(selections)
    for path, receipt_key, root in [
        ('metadata/selected_taxon_inputs_readback_20260927.json', 'producer_receipt_sha256', source),
        ('metadata/species_distance_kernel_readback_20260927.json', 'source_receipt_sha256', kernels)]:
        proof = json.loads(Path(path).read_text())
        assert proof['status'].startswith('passed')
        assert proof[receipt_key] == sha(root / 'receipt.json')
        pins[path] = sha(path)
    nodes = pd.read_csv(source / 'selected_node_taxa.tsv', sep='\t')
    target = nodes[nodes.role.eq('target')].set_index('node_id')
    background = nodes[nodes.role.eq('background')].set_index('node_id')
    selected = pd.read_csv(selections, sep='\t', usecols=['target_id', 'background_id'])
    pairs = selected.groupby(['target_id', 'background_id']).size().rename('selected_record_uses').reset_index()
    assert pairs.selected_record_uses.sum() == len(selected) == 2786912
    for key in ['guide', 'family']:
        assert pairs.target_id.map(target[key]).eq(pairs.background_id.map(background[key])).all()
    pairs['focal_index'] = pairs.target_id.map(target.kernel_index_a)
    pairs['background_index_a'] = pairs.background_id.map(background.kernel_index_a)
    pairs['background_index_b'] = pairs.background_id.map(background.kernel_index_b)
    taxa = json.loads((kernels / 'taxa.json').read_text())
    patterns = {}
    assignments = []
    for a, b, c in pairs[['focal_index', 'background_index_a', 'background_index_b']].itertuples(index=False, name=None):
        weights = {}
        for index, weight in [(a, 2), (b, -1), (c, -1)]:
            weights[int(index)] = weights.get(int(index), 0) + weight
        pattern = tuple(sorted((i, w) for i, w in weights.items() if w))
        assert sum(w for i, w in pattern) == 0
        text = json.dumps(pattern, separators=(',', ':'))
        pattern_id = hashlib.sha256(text.encode()).hexdigest()
        assert pattern_id not in patterns or patterns[pattern_id] == pattern
        patterns[pattern_id] = pattern
        assignments.append(pattern_id)
    pairs['species_pattern_id'] = assignments
    labels = sorted(patterns)
    rows, cols, data = [], [], []
    pattern_rows = []
    for r, label in enumerate(labels):
        pattern = patterns[label]
        pattern_rows.append(dict(species_pattern_id=label, row_index=r,
                                 taxa=json.dumps([taxa[i] for i,w in pattern], separators=(',', ':')),
                                 twice_weights=json.dumps([w for i,w in pattern], separators=(',', ':'))))
        for i, w in pattern:
            rows.append(r); cols.append(i); data.append(w / 2)
    design = sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(taxa)))
    assert np.all(np.asarray(design.sum(axis=1)) == 0)
    out.mkdir(exist_ok=False)
    sparse.save_npz(out / 'species_contrast_design.npz', design)
    pd.DataFrame(pattern_rows).to_csv(out / 'patterns.tsv', sep='\t', index=False)
    pairs.to_csv(out / 'selected_pair_patterns.tsv.gz', sep='\t', index=False)
    (out / 'taxa.json').write_text(json.dumps(taxa, indent=2) + '\n')
    # Every exported pattern is checked through independent pair-distance algebra.
    variance_rows = []
    pattern_index = {label: i for i,label in enumerate(labels)}
    pair_rows = pairs.species_pattern_id.map(pattern_index).to_numpy()
    a,b,c = [pairs[k].to_numpy(dtype=int) for k in ['focal_index', 'background_index_a', 'background_index_b']]
    maximum_error = 0.
    for path in sorted(kernels.glob('*.npz')):
        with np.load(path) as arrays:
            covariance, distance = arrays['centered_kernel'], arrays['distances']
            diagonal = np.asarray(design.multiply(design @ covariance).sum(axis=1)).ravel()
            independent = .5 * distance[a,b] + .5 * distance[a,c] - .25 * distance[b,c]
        np.testing.assert_allclose(diagonal[pair_rows], independent, atol=1e-12, rtol=1e-10)
        assert diagonal.min() >= -1e-12
        maximum_error = max(maximum_error, float(np.max(abs(diagonal[pair_rows] - independent))))
        variance_rows.extend(dict(tree=path.stem, species_pattern_id=label, kernel_quadratic_form=float(value)) for label,value in zip(labels,diagonal))
    pd.DataFrame(variance_rows).to_csv(out / 'kernel_quadratic_forms.tsv', sep='\t', index=False)
    # Independent raw-record sparse reconstruction, in bounded chunks.
    reconstructed = sparse.load_npz(out / 'species_contrast_design.npz')
    assert reconstructed.shape == design.shape
    exported = pd.read_csv(out / 'selected_pair_patterns.tsv.gz', sep='\t')
    pd.testing.assert_frame_equal(exported, pairs)
    for start in range(0, len(pairs), 10000):
        stop = min(start + 10000, len(pairs)); n = stop-start
        expected = sparse.coo_matrix((np.tile([1., -.5, -.5], n),
                                     (np.repeat(np.arange(n), 3), np.column_stack([a[start:stop],b[start:stop],c[start:stop]]).ravel())), shape=(n,526)).tocsr()
        delta = reconstructed[pair_rows[start:stop]] - expected
        assert delta.nnz == 0
    for path,digest in pins.items():
        assert sha(path) == digest, path
    result = dict(status='complete_matched_species_contrast_design_with_full_pair_readback',
                  selected_records=len(selected), unique_selected_pairs=len(pairs), patterns=len(patterns),
                  zero_patterns=sum(not p for p in patterns.values()), trees=5,
                  pair_tree_quadratic_forms_checked=len(pairs)*5, maximum_quadratic_form_error=maximum_error,
                  pins=pins, script_sha256=sha(__file__), artifacts={p.name:sha(p) for p in out.iterdir()},
                  scope='Working additive endpoint contrast W: focal species minus the mean of background endpoint species, with exact cancellation for shared taxa. Retains every selected pair and source record multiplicity before structural qualification. W C W-transpose is a prospective nuisance covariance, not an estimated structural process, effect or dated rate. Quadratic forms are unscaled kernel quantities, not fitted variances. No model fitting or GPU inference.')
    (out / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts']}, indent=2))


if __name__ == '__main__':
    main()
