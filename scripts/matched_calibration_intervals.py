"""Fixed-sample Monte Carlo coverage intervals with unresolved outcomes retained."""
import math
from scipy.stats import beta
from refit_matched_simulation import coverage_accounting


def exact_coverage_bounds(covered, unresolved, attempted, alpha=.05):
    """Envelope of exact binomial intervals over every unresolved completion.

    Requires a prespecified number of independent replicates under the same
    generating model. These are marginal Monte Carlo intervals, not biological
    coefficient intervals or intervals valid under optional stopping.
    """
    if any(type(n) is not int for n in [covered, unresolved, attempted]):
        raise ValueError('Integer outcome counts required')
    if attempted < 1 or covered < 0 or unresolved < 0 or covered+unresolved > attempted:
        raise ValueError('Invalid attempted/covered/unresolved counts')
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('Alpha must lie strictly between zero and one')
    maximum_covered = covered+unresolved
    lower = 0. if covered == 0 else float(beta.ppf(alpha/2, covered, attempted-covered+1))
    upper = 1. if maximum_covered == attempted else float(beta.isf(alpha/2, maximum_covered+1, attempted-maximum_covered))
    assert 0 <= lower <= upper <= 1
    return dict(attempted=attempted, covered=covered, unresolved=unresolved,
                alpha=alpha, lower=lower, upper=upper,
                observed_fraction_lower=covered/attempted,
                observed_fraction_upper=maximum_covered/attempted)


def summarize_replicates(records, coefficients, alpha=.05):
    accounting = coverage_accounting(records, coefficients)
    return dict(accounting=accounting, marginal_monte_carlo_intervals=[
        exact_coverage_bounds(covered, accounting['unresolved'], accounting['attempted'], alpha)
        for covered in accounting['covered_qualified']])
