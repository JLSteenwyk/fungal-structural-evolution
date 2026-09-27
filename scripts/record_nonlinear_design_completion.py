"""Close full nonlinear design verification and summarize numerical conditioning."""
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    states = {}
    for unit in ['fungal-nonlinear-matched-designs-20260927.service',
                 'fungal-nonlinear-matched-design-readback-20260927.service']:
        state = dict(line.split('=', 1) for line in subprocess.check_output(
            ['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '-p', 'Result', '-p', 'ExecMainStatus'], text=True).splitlines())
        assert state == dict(ActiveState='inactive', Result='success', ExecMainStatus='0'), (unit, state)
        states[unit] = state
    pp = Path('metadata/nonlinear_matched_design_plan_20260927.json')
    ap = Path('metadata/nonlinear_matched_design_readback_20260927.json')
    plan, audit = json.loads(pp.read_text()), json.loads(ap.read_text())
    root = Path(plan['output'])
    rp = root/'receipt.json'
    receipt = json.loads(rp.read_text())
    assert receipt['plan_sha256'] == sha(pp)
    assert audit['status'] == 'passed_full_nonlinear_matched_design_diagnostics_readback'
    assert audit['source_receipt_sha256'] == sha(rp)
    for path, digest in plan['pins'].items():
        assert sha(path) == digest, path
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest, name
    table = pd.read_csv(root/'design_diagnostics.tsv', sep='\t')
    assert len(table) == audit['strata'] == receipt['strata'] == 165888
    keys = ['guide', 'policy', 'scenario_id', 'boundary', 'mask', 'cohort', 'screen', 'target_order', 'background_order']
    assert not table.duplicated(keys+['polynomial_degree']).any()
    assert set(table.polynomial_degree) == {2, 3}
    summaries = []
    for degree, frame in table.groupby('polynomial_degree'):
        assert len(frame) == 82944
        condition = frame.scaled_condition.to_numpy(float)
        spectra = [np.array(json.loads(s)) for s in frame.scaled_singular_values]
        for row, singular in zip(frame.itertuples(index=False), spectra):
            rank = int((singular > receipt['rank_absolute_tolerance']).sum())
            assert rank == row.centered_rank
            assert row.full_rank_after_constant_removal == int(rank == row.active_covariates)
        finite = condition[np.isfinite(condition)]
        assert len(finite) == int(frame.full_rank_after_constant_removal.sum())
        quantiles = np.quantile(finite, [0, .5, .95, 1])
        summaries.append(dict(polynomial_degree=int(degree), settings=len(frame),
            rank_deficient=int(frame.full_rank_after_constant_removal.eq(0).sum()),
            scaled_condition_min=float(quantiles[0]), scaled_condition_median=float(quantiles[1]),
            scaled_condition_p95=float(quantiles[2]), scaled_condition_max=float(quantiles[3]),
            minimum_scaled_singular_value=min(float(s[-1]) for s in spectra),
            outside_some_marginal_range=int(frame.zero_in_each_marginal_range.eq(0).sum())))
    paired = table.pivot(index=keys, columns='polynomial_degree', values='scaled_condition')
    assert len(paired) == 82944 and not paired.isna().any().any()
    # Adding a standardized column cannot improve the spectral condition number.
    assert (paired[3] >= paired[2]-1e-8).all()
    ratio = paired[3]/paired[2]
    result = dict(status='complete_verified_nonlinear_design_diagnostics',
                  terminal_states=states, source_receipt_sha256=sha(rp), audit_sha256=sha(ap),
                  plan_sha256=sha(pp), script_sha256=sha(__file__), summaries=summaries,
                  cubic_to_quadratic_condition_ratio_median=float(ratio.median()),
                  cubic_to_quadratic_condition_ratio_maximum=float(ratio.max()),
                  scope='All 165888 design rows passed independent SQL-moment readback. Condition numbers describe centered population-SD-scaled covariates before covariance weighting; repeated sensitivity settings are not independent samples. Full rank does not establish effect precision, nonlinear reference support, fitted model adequacy or calibrated inference.')
    Path('metadata/nonlinear_matched_design_completed_20260927.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
