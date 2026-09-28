"""Generate Gaussian working-model responses without forming an n-by-n covariance.

Simulation is conditional on a fixed design, species factor and variance
components. Calibration still requires refitting every simulated response.
"""
import numpy as np


def simulate(background, family, species_factor, design, beta, scale, ratios,
             draws, rng):
    background, family = np.asarray(background), np.asarray(family)
    factor, x, beta = (np.asarray(v, dtype=float) for v in (species_factor, design, beta))
    ratios = np.asarray(ratios, dtype=float)
    if background.ndim != 1 or family.shape != background.shape or not len(background):
        raise ValueError('Expected equally sized nonempty group label vectors')
    n = len(background)
    if factor.ndim != 2 or factor.shape[0] != n or not np.isfinite(factor).all():
        raise ValueError('Invalid species factor')
    if x.ndim != 2 or x.shape[0] != n or not x.shape[1] or beta.shape != (x.shape[1],):
        raise ValueError('Invalid fixed-effect dimensions')
    if not np.isfinite(x).all() or not np.isfinite(beta).all():
        raise ValueError('Nonfinite fixed effects')
    if not np.isfinite(scale) or scale <= 0 or ratios.shape != (3,) or not np.isfinite(ratios).all() or (ratios < 0).any():
        raise ValueError('Positive scale and three nonnegative variance ratios required')
    if type(draws) is not int or draws < 1:
        raise ValueError('Positive integer draw count required')
    _, bg = np.unique(background, return_inverse=True)
    _, fam = np.unique(family, return_inverse=True)
    ng, nf = int(bg.max())+1, int(fam.max())+1
    membership = {}
    for b, f in zip(bg, fam):
        if b in membership and membership[b] != f:
            raise ValueError('Background groups must be nested within family components')
        membership[b] = f
    # Keep draw ordering fixed even for zero variance components. Generator
    # objects are supplied by callers so seeds and independent streams are explicit.
    residual = rng.standard_normal((n, draws))
    bg_effect = rng.standard_normal((ng, draws))
    family_effect = rng.standard_normal((nf, draws))
    species_effect = rng.standard_normal((factor.shape[1], draws))
    noise = residual + np.sqrt(ratios[0])*bg_effect[bg] + np.sqrt(ratios[1])*family_effect[fam]
    if factor.shape[1]:
        noise += np.sqrt(ratios[2])*(factor @ species_effect)
    return (x @ beta)[:, None] + np.sqrt(scale)*noise
