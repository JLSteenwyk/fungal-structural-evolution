"""Check that equal-likelihood parameter variation survives summarization."""
import json
from pathlib import Path
from summarize_mg94_multistart_parameters import parameter_ranges
from readback_whole_proteome_catalog import sha
rows=[{'log_likelihood':-10.},{'log_likelihood':-10.000001},{'log_likelihood':-11.}]
pars=[{'omega':.2,'branch':0.},{'omega':.8,'branch':2.},{'omega':5.,'branch':100.}]
r=parameter_ranges(rows,pars,1e-5)
assert r[0]['near_best_max']==2 and r[0]['all_start_max']==100
assert r[1]['near_best_min']==.2 and r[1]['near_best_max']==.8 and r[1]['all_start_max']==5
assert all(x['near_best_starts']==2 for x in r)
r=parameter_ranges(rows,pars,0)
assert all(x['near_best_starts']==1 and x['near_best_range']==0 for x in r)
try:
    parameter_ranges(rows,[{'omega':float('nan'),'branch':0.},*pars[1:]],1e-5)
except AssertionError:
    pass
else:
    raise AssertionError('Nonfinite parameter accepted')
result={'status':'passed_tied_likelihood_parameter_variation_and_poor_start_exclusion_checks',
        'script_sha256':sha('scripts/summarize_mg94_multistart_parameters.py'),
        'checker_sha256':sha(__file__),
        'checks':['near-best parameter variability retained','poorer likelihood start excluded only from near-best range','all-start extrema preserved','single optimum produces zero descriptive range','nonfinite parameter rejected']}
Path('metadata/mg94_parameter_range_checks_20260928.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'])
