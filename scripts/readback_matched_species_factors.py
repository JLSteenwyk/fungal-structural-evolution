#!/usr/bin/env python3
"""Check all factor covariances against independently parsed tree edge paths."""
import json
from pathlib import Path
import dendropy
import numpy as np
from scipy import sparse
from screen_duplication_domain_alignment_coverage import sha


def main():
    root = Path('results/phylogeny/matched-species-covariance-factors-20260927-v1')
    receipt = json.loads((root / 'receipt.json').read_text())
    for path,digest in receipt['pins'].items():
        assert sha(path) == digest, path
    for name,digest in receipt['artifacts'].items():
        assert sha(root / name) == digest, name
    plan_path = Path('metadata/species_distance_kernel_plan_20260927.json')
    plan = json.loads(plan_path.read_text())
    kernel_receipt = json.loads(Path('results/phylogeny/species-distance-kernels-20260927-v1/receipt.json').read_text())
    assert sha(plan_path) == kernel_receipt['plan_sha256']
    taxa = json.loads((root / 'taxa.json').read_text())
    index = {label:i for i,label in enumerate(taxa)}
    design_root = Path('results/phylogeny/matched-species-contrasts-20260927-v1')
    assert (root / 'patterns.tsv').read_bytes() == (design_root / 'patterns.tsv').read_bytes()
    assert (root / 'taxa.json').read_bytes() == (design_root / 'taxa.json').read_bytes()
    design = sparse.load_npz(design_root / 'species_contrast_design.npz')
    q = np.load(root / 'pattern_basis.npy')
    projection = np.load(root / 'species_projection.npy')
    assert q.shape == (4568,242) and projection.shape == (242,526)
    np.testing.assert_allclose(q.T@q, np.eye(242), atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(q@projection, design.toarray(), atol=1e-12, rtol=1e-10)
    summaries = {}
    for label,spec in plan['trees'].items():
        assert sha(spec['tree']) == kernel_receipt['trees'][label]['tree_sha256']
        tree = dendropy.Tree.get(path=spec['tree'], schema='newick', preserve_underscores=True)
        assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == set(taxa)
        edges = [node for node in tree.preorder_node_iter() if node.parent_node is not None]
        features = np.zeros((len(taxa), len(edges)))
        for column,node in enumerate(edges):
            length = float(node.edge_length)
            assert np.isfinite(length) and length >= 0
            for leaf in node.leaf_iter():
                features[index[leaf.taxon.label], column] = np.sqrt(length)
        paths = design @ features
        with np.load(root / (label+'.npz')) as arrays:
            factor,core,lower = arrays['factor'],arrays['reduced_covariance'],arrays['cholesky']
        assert factor.shape == (4568,242)
        np.testing.assert_array_equal(np.triu(lower,1), 0)
        assert np.all(np.diag(lower) > 0)
        np.testing.assert_allclose(factor, q@lower, atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(lower@lower.T, core, atol=1e-11, rtol=1e-10)
        projected_paths = projection @ features
        np.testing.assert_allclose(projected_paths@projected_paths.T, core, atol=1e-10, rtol=1e-10)
        maximum = 0.
        count = 0
        for start in range(0,len(factor),128):
            end = min(start+128,len(factor))
            expected = paths[start:end]@paths.T
            observed = factor[start:end]@factor.T
            np.testing.assert_allclose(observed, expected, atol=1e-11, rtol=1e-10)
            maximum = max(maximum,float(np.max(abs(observed-expected))))
            count += expected.size
        assert count == receipt['trees'][label]['covariance_entries_checked']
        np.testing.assert_allclose(np.linalg.eigvalsh(core)[0], receipt['trees'][label]['minimum_core_eigenvalue'], atol=1e-12, rtol=1e-8)
        summaries[label] = dict(covariance_entries_checked=count, maximum_tree_edge_error=maximum)
        print('Checked full covariance from tree edges:', label, flush=True)
    for path,digest in receipt['pins'].items():
        assert sha(path) == digest, path
    result = dict(status='passed_full_matched_species_factor_tree_readback',
                  trees=summaries, covariance_entries_checked=sum(x['covariance_entries_checked'] for x in summaries.values()),
                  source_receipt_sha256=sha(root / 'receipt.json'), script_sha256=sha(__file__),
                  scope='Every ordered pattern-pair covariance reconstructed from independently parsed DendroPy tree edge paths. Complete basis, species projection, Cholesky and core covariance checked. No fitted variance multiplier, inverse, regularization or evolutionary effect.')
    Path('metadata/matched_species_factor_readback_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()
