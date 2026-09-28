"""Audit every saved paired information result after successful runner exit."""
import json
from pathlib import Path
import numpy as np
from ancestral_chain_attempt import sha, write_json


def main():
    root = Path('results/model_validation/matched-information-shortcut-comparison-20260928-v1')
    source = root / 'receipt.json'
    r = json.loads(source.read_text())
    assert r['status'] == 'passed_all_15_real_input_information_shortcut_comparisons'
    for path, digest in r['source_pins'].items():
        assert sha(path) == digest, path
    for name, digest in r['artifacts'].items():
        assert sha(root / name) == digest, name
    reference_root = Path('results/model_validation/matched-information-timing-20260928-v1')
    ref = json.loads((reference_root / 'receipt.json').read_text())
    assert len(r['jobs']) == len(ref['jobs']) == 15
    assert [j['job'] for j in r['jobs']] == [j['job'] for j in ref['jobs']]
    maximum = 0.
    for fast, original in zip(r['jobs'], ref['jobs']):
        assert fast['reference_seconds'] == original['elapsed_seconds']
        assert fast['shortcut_seconds'] > 0 and np.isfinite(fast['shortcut_seconds'])
        with np.load(root / fast['artifact'], allow_pickle=False) as a, np.load(reference_root / original['artifact'], allow_pickle=False) as b:
            assert set(a.files) == set(b.files)
            for key in a.files:
                np.testing.assert_allclose(a[key], b[key], rtol=1e-7, atol=1e-8)
                if a[key].dtype.kind != 'b':
                    error = float(np.max(abs(a[key]-b[key])))
                    assert error == fast['maximum_absolute_errors'][key]
            maximum = max(maximum, fast['maximum_absolute_errors']['expected_reml_information'])
            info = a['expected_reml_information']
            diagonal = np.sqrt(np.diag(info))
            eig = np.linalg.eigvalsh(info / np.outer(diagonal, diagonal))
            np.testing.assert_allclose(eig, a['normalized_information_eigenvalues'], rtol=1e-10, atol=1e-12)
            tolerance = 4*np.finfo(float).eps*max(float(np.max(abs(eig))), 1.)
            assert int(np.sum(eig > tolerance)) == fast['numerical_information_rank']
    before = sum(j['reference_seconds'] for j in r['jobs'])
    after = sum(j['shortcut_seconds'] for j in r['jobs'])
    assert before == r['reference_seconds_sum'] and after == r['shortcut_seconds_sum']
    result = dict(status='passed_all_saved_information_comparison_artifacts', cases=15,
        source_receipt_sha256=sha(source), script_sha256=sha(__file__),
        maximum_information_absolute_disagreement=maximum,
        reference_seconds=before, shortcut_seconds=after,
        measured_time_reduction_fraction=1-after/before,
        ranks=[j['numerical_information_rank'] for j in r['jobs']],
        scope='All saved comparison arrays and timing arithmetic checked. Single timings; '
              'selected cases do not establish global speed, coverage, or boundary regularity.')
    write_json(Path('metadata/matched_information_shortcut_readback_20260928.json'), result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
