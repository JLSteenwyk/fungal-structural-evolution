"""Check descriptive weighting, missing dispositions and invalid input handling."""
import numpy as np
import pandas as pd
from summarize_selected_matched_residuals import summarize, GOOD

rows = []
for identifier, n, a, b in [('a', 10, 4, 1), ('b', 100, 2, 0)]:
    rows.append(dict(fit_input_id=identifier, tree='one', status=GOOD, records=n,
                     absolute_above_2=a, absolute_above_3=b, mean=0.,
                     second_raw_moment=1., third_raw_moment=0., fourth_raw_moment=3.,
                     maximum_fixed_effect_variance_fraction=.1))
rows.append(dict(fit_input_id='c', tree='one', status='residual_computation_requires_review'))
frame = pd.DataFrame(rows)
out, counts, summary = summarize(frame)
assert len(out) == 3 and counts.fits.sum() == 3 and counts.tree_total_fits.eq(3).all()
record = summary[summary.metric.eq('fraction_absolute_above_2')].iloc[0]
assert record.available_fits == 2 and record.unavailable_fits == 1
assert np.isclose(record['median'], .21) and not np.isclose(record['median'], 6/110)
assert np.isnan(out.iloc[2].fraction_absolute_above_2)
for variant in ['duplicate', 'impossible_tail', 'nonfinite']:
    bad = frame.copy()
    if variant == 'duplicate':
        bad = pd.concat([bad, bad.iloc[[0]]], ignore_index=True)
    if variant == 'impossible_tail':
        bad.loc[0, 'absolute_above_3'] = 5
    if variant == 'nonfinite':
        bad.loc[0, 'fourth_raw_moment'] = np.inf
    try:
        summarize(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError(variant)
print('Passed weighting, denominators, retained failures and malformed-input checks')
