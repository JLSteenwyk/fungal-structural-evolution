#!/usr/bin/env python3
"""Reconstruct exact inputs for a frozen optimization-review snapshot."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from screen_duplication_alignment_reuse import sha


def main():
    pp = Path('metadata/full_polynomial_ml_plan_20260927.json')
    plan = json.loads(pp.read_text())
    census = Path('results/model_validation/polynomial-optimization-flag-census-20260928-v1')
    cr = json.loads((census / 'receipt.json').read_text())
    bindings = {str(pp): cr['source_plan_sha256'], str(census / 'receipt.json'): sha(census / 'receipt.json')}
    bindings.update(cr['source_hashes'])
    bindings.update({str(census / n): h for n, h in cr['artifacts'].items()})
    flags = pd.read_csv(census / 'flagged_fits.tsv', sep='\t')
    ids = set(flags.fit_input_id)
    recipes = {}
    for key in ['inventory', 'nonlinear_inventory']:
        folder = Path(plan[key])
        receipt = json.loads((folder / 'receipt.json').read_text())
        audit_path = plan[key + '_audit']
        audit = json.loads(Path(audit_path).read_text())
        assert audit['source_receipt_sha256'] == sha(folder / 'receipt.json') and audit['status'].startswith('passed_full_')
        bindings[audit_path] = sha(audit_path)
        bindings[str(folder / 'receipt.json')] = sha(folder / 'receipt.json')
        rp = folder / 'unique_fit_recipes.jsonl'
        bindings[str(rp)] = receipt['artifacts'][rp.name]
        for row in map(json.loads, rp.open()):
            if row['fit_input_id'] in ids:
                assert row['fit_input_id'] not in recipes
                recipes[row['fit_input_id']] = dict(row, polynomial_degree=row.get('polynomial_degree', 1))
    assert set(recipes) == ids and len(ids) == 449 and len(flags) == 601
    for key in ['nodes', 'pairs', 'selections']:
        bindings[plan[key]] = plan['pins'][plan[key]]
    source = Path(plan['summaries'])
    nonlinear = Path(plan['nonlinear'])
    for folder in [source, nonlinear, Path(plan['factors'])]:
        receipt = json.loads((folder / 'receipt.json').read_text())
        bindings[str(folder / 'receipt.json')] = plan['pins'].get(str(folder / 'receipt.json'), sha(folder / 'receipt.json'))
        bindings.update({str(folder / n): h for n, h in receipt['artifacts'].items()})
    def verify():
        for path, digest in bindings.items():
            assert sha(path) == digest, path
    verify()
    nodes = pd.read_csv(plan['nodes'], sep='\t')
    targets = nodes.loc[nodes.role.eq('target'), ['node_id', 'guide', 'family_component']].rename(columns={'node_id': 'target_id'})
    pairs = pd.read_csv(plan['pairs'], sep='\t', usecols=['target_id', 'background_id', 'species_pattern_id']).merge(targets, on='target_id', validate='many_to_one')
    pairs['row_identity'] = [hashlib.sha256(json.dumps(list(row), separators=(',', ':')).encode()).hexdigest() for row in pairs[['target_id', 'background_id', 'family_component', 'species_pattern_id']].itertuples(index=False, name=None)]
    selected = pd.read_csv(plan['selections'], sep='\t', usecols=['target_id', 'background_id', 'policy', 'scenario_id', 'domain_config_id']).merge(pairs, on=['target_id', 'background_id'], validate='many_to_one').sort_values('target_id', kind='stable')
    assert len(selected) == 2786912
    patterns = pd.read_csv(Path(plan['factors']) / 'patterns.tsv', sep='\t').set_index('species_pattern_id').row_index
    parts = json.loads((source / 'partition_manifest.json').read_text())
    numeric = ['rmsd_difference', 'identity_difference', 'original_coverage_difference', 'log_aligned_length_ratio', 'confidence_fraction_difference', 'identity_power_2_difference', 'identity_power_3_difference']
    out = Path('results/model_validation/flagged-polynomial-inputs-20260928-v1')
    out.mkdir(exist_ok=False)
    exported = []
    for part in parts:
        needed = [r for r in recipes.values() if r['partition_index'] == part['index']]
        if not needed:
            continue
        assert sha(source / part['path']) == part['sha256']
        values = pd.read_parquet(source / part['path'], columns=['domain_config_id'] + numeric[:5])
        extra = pd.read_parquet(nonlinear / f"{part['index']:03d}.parquet")
        assert set(values.domain_config_id) == set(extra.domain_config_id)
        values = values.merge(extra[['domain_config_id'] + numeric[5:]], on='domain_config_id', validate='one_to_one')
        records = selected.merge(values, on='domain_config_id', validate='many_to_one', sort=False)
        groups = records.groupby(['guide', 'policy', 'scenario_id']).indices
        for recipe in needed:
            frame = records.iloc[groups[tuple(recipe[k] for k in ['guide', 'policy', 'scenario_id'])]]
            assert frame.target_id.is_monotonic_increasing and len(frame) == recipe['records']
            matrix = np.ascontiguousarray(frame[numeric[:recipe['polynomial_degree'] + 4]].to_numpy(), dtype='<f8')
            matrix[matrix == 0] = 0.
            identities = np.asarray(frame.row_identity, dtype='S64')
            assert hashlib.sha256(matrix.tobytes()).hexdigest() == recipe['values_sha256']
            assert hashlib.sha256(identities.tobytes()).hexdigest() == recipe['ordered_identity_sha256']
            mapped = frame.species_pattern_id.map(patterns)
            assert mapped.notna().all() and np.isfinite(matrix).all()
            arrays = dict(matrix=matrix, row_identity=identities,
                          background=pd.factorize(frame.background_id, sort=True)[0],
                          family=pd.factorize(frame.family_component, sort=True)[0],
                          pattern_rows=mapped.to_numpy(dtype=int))
            name = recipe['fit_input_id'] + '.npz'
            np.savez_compressed(out / name, **arrays)
            with np.load(out / name) as saved:
                assert set(saved.files) == set(arrays)
                for k, value in arrays.items():
                    np.testing.assert_array_equal(saved[k], value)
            exported.append(dict(recipe=recipe, path=name, sha256=sha(out / name)))
        print('Reconstructed flagged partition', part['index'], 'inputs', len(exported), '/449', flush=True)
    assert len(exported) == 449
    (out / 'input_manifest.json').write_text(json.dumps(exported, indent=2) + '\n')
    flags.to_csv(out / 'flagged_tree_fits.tsv', sep='\t', index=False)
    verify()
    receipt = dict(status='complete_flagged_polynomial_input_reconstruction_pending_independent_readback', inputs=len(exported), tree_fits=len(flags), record_occurrences=sum(r['recipe']['records'] for r in exported), source_hashes=bindings, script_sha256=sha(__file__), artifacts={p.name: sha(p) for p in out.iterdir()}, scope='Frozen 601-flag snapshot only. All 449 exact numeric and ordered-identity recipe hashes matched; all exported arrays read back. No production fits modified or flags resolved; independent input validation and fit resource estimate required.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ['source_hashes', 'artifacts']}), flush=True)


if __name__ == '__main__':
    main()
