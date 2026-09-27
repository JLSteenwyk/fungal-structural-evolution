"""Link fully audited polynomial inputs to identical-observation linear references."""
import json
from pathlib import Path
import subprocess
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    units = ['fungal-nonlinear-model-input-inventory-20260927.service',
             'fungal-nonlinear-model-input-readback-20260927.service']
    states = {}
    for unit in units:
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), (unit, state)
        states[unit] = state
    plans = [Path('metadata/full_matched_fit_inventory_plan_20260927.json'),
             Path('metadata/nonlinear_model_input_inventory_plan_20260927.json')]
    audits = [Path('metadata/full_matched_fit_inventory_readback_20260927.json'),
              Path('metadata/nonlinear_model_input_inventory_readback_20260927.json')]
    statuses = ['passed_full_matched_fit_inventory_readback', 'passed_full_nonlinear_model_input_inventory_readback']
    tables, recipes, configs, pins = [], [], [], {}
    for pp, ap, status in zip(plans, audits, statuses):
        config = json.loads(pp.read_text())
        root = Path(config['output'])
        rp = root/'receipt.json'
        receipt, audit = json.loads(rp.read_text()), json.loads(ap.read_text())
        assert audit['status'] == status and audit['source_receipt_sha256'] == sha(rp)
        assert receipt['plan_sha256'] == sha(pp)
        for name, digest in receipt['artifacts'].items():
            assert sha(root/name) == digest
            pins[str(root/name)] = digest
        pins.update({str(pp):sha(pp), str(ap):sha(ap), str(rp):sha(rp)})
        tables.append(pd.read_csv(root/'full_setting_fit_map.tsv', sep='\t'))
        entries = [json.loads(line) for line in (root/'unique_fit_recipes.jsonl').open()]
        lookup = {row['fit_input_id']:row for row in entries}
        assert len(lookup) == len(entries) == receipt['unique_record_inputs']
        recipes.append(lookup)
        configs.append(config)
    for key in ['summaries', 'nodes', 'pairs', 'selections']:
        assert configs[0][key] == configs[1][key]
    for path in [configs[0][k] for k in ['nodes', 'pairs', 'selections']] + [configs[0]['summaries']+'/receipt.json']:
        assert configs[0]['pins'][path] == configs[1]['pins'][path] == sha(path)
        pins[path] = sha(path)
    keys = ['partition_index', 'guide', 'policy', 'scenario_id', 'boundary', 'mask', 'cohort', 'screen', 'target_order', 'background_order', 'records']
    linear, nonlinear = tables
    assert len(linear) == 82944 and len(nonlinear) == 165888
    assert not linear.duplicated(keys).any() and not nonlinear.duplicated(keys+['polynomial_degree']).any()
    wide = linear.rename(columns={'fit_input_id':'linear_input_id'})
    for degree, name in [(2, 'quadratic_input_id'), (3, 'cubic_input_id')]:
        candidate = nonlinear[nonlinear.polynomial_degree.eq(degree)].drop(columns='polynomial_degree').rename(columns={'fit_input_id':name})
        wide = wide.merge(candidate, on=keys, how='outer', validate='one_to_one', indicator=True)
        assert wide['_merge'].eq('both').all()
        wide = wide.drop(columns='_merge')
    assert len(wide) == 82944
    columns = ['linear_input_id', 'quadratic_input_id', 'cubic_input_id']
    unique = wide[columns].drop_duplicates().sort_values(columns).reset_index(drop=True)
    for row in unique.itertuples(index=False):
        base = recipes[0][row.linear_input_id]
        for degree, identifier in [(2, row.quadratic_input_id), (3, row.cubic_input_id)]:
            expanded = recipes[1][identifier]
            assert expanded['polynomial_degree'] == degree
            assert expanded['records'] == base['records']
            assert expanded['ordered_identity_sha256'] == base['ordered_identity_sha256']
            assert expanded['numeric_columns'][:len(base['numeric_columns'])] == base['numeric_columns']
    out = Path('results/model_validation/polynomial-model-input-links-20260927-v1')
    out.mkdir(parents=True, exist_ok=False)
    artifacts = {}
    for name, table in [('full_setting_links.tsv', wide), ('unique_comparison_sets.tsv', unique)]:
        path = out/name
        table.to_csv(path, sep='\t', index=False)
        pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'), table.reset_index(drop=True))
        artifacts[name] = sha(path)
    for path, digest in pins.items():
        assert sha(path) == digest
    result = dict(status='complete_verified_identical_observation_polynomial_input_links',
                  terminal_states=states, settings=len(wide), unique_comparison_sets=len(unique),
                  linear_inputs_with_multiple_expanded_sets=int(unique.groupby('linear_input_id').size().gt(1).sum()),
                  source_hashes=pins, artifacts=artifacts, script_sha256=sha(__file__),
                  scope='All settings joined one-to-one with exact record counts and ordered observation identities; audited shared original data establish common response/base covariates. Input links only. Existing REML objectives are not ML reference fits, and no likelihood comparison or significance test is performed.')
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    Path('metadata/polynomial_model_input_links_completed_20260927.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','settings','unique_comparison_sets','linear_inputs_with_multiple_expanded_sets']}))


if __name__ == '__main__':
    main()
