"""Known hull cases and deliberate certificate corruption checks."""
import copy
import json
from pathlib import Path
from matched_joint_covariate_support import assess_support
from check_joint_support_certificate import check_certificate
from screen_duplication_domain_alignment_coverage import sha

cases = [([[-1,-1],[-1,1],[1,-1],[1,1]], 'zero_supported_to_numeric_tolerance'),
         ([[0,0],[1,0],[0,1]], 'zero_supported_to_numeric_tolerance'),
         ([[1,-.1],[-.1,1],[1,1]], 'zero_outside_joint_convex_hull'),
         ([[1,0],[1,1]], 'zero_outside_joint_convex_hull'),
         ([[0,0],[0,0]], 'zero_supported_to_numeric_tolerance')]
for x, expected in cases:
    result = assess_support(x)
    assert result['classification'] == expected
    check_certificate(x, json.loads(json.dumps(result)))
x = cases[2][0]
result = assess_support(x)
changes = {
    'scales': [2,2], 'support_weights': [0]*len(result['support_weights']),
    'support_indices': [999]*len(result['support_indices']),
    'separating_direction': [0,0], 'classification': 'zero_supported_to_numeric_tolerance',
    'barycenter': [99,99], 'certificate_valid': False,
    'primal_distance': 99., 'solver_objective': 99.,
}
for key, value in changes.items():
    bad = copy.deepcopy(result)
    bad[key] = value
    try:
        check_certificate(x, bad)
    except (AssertionError, ValueError):
        continue
    raise AssertionError('Accepted corruption: '+key)
path = Path('metadata/joint_support_certificate_checks_20260927.json')
path.write_text(json.dumps(dict(status='passed', known_cases=len(cases), rejected_corruptions=len(changes),
    pins={p:sha(p) for p in [__file__,'scripts/check_joint_support_certificate.py','scripts/matched_joint_covariate_support.py']}),indent=2)+'\n')
print(path.read_text())
