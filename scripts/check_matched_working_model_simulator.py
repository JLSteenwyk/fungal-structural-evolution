"""Exact finite latent-basis covariance checks for future calibration simulation."""
import itertools
from pathlib import Path
import numpy as np
from ancestral_chain_attempt import sha, write_json
from simulate_matched_working_model import simulate


class BasisDraws:
    def __init__(self, dimensions):
        self.values = np.sqrt(dimensions)*np.concatenate([np.eye(dimensions), -np.eye(dimensions)], axis=1)
        self.offset = 0

    def standard_normal(self, shape):
        rows, draws = shape
        assert draws == self.values.shape[1]
        value = self.values[self.offset:self.offset+rows]
        assert value.shape == shape
        self.offset += rows
        return value


def main():
    bg = np.array(['a','a','b','c','c','d'])
    fam = np.array(['x','x','x','y','y','y'])
    x = np.column_stack([np.ones(6), np.arange(6)/5])
    beta = np.array([.7, -.3]); scale = 2.4
    factor = np.array([[1,0],[.8,.2],[.5,.3],[0,1],[.2,.9],[.1,.5]])
    maximum = 0.; cases = 0
    for rank in [0, 2]:
        f = factor[:, :rank]
        dimensions = len(bg)+len(set(bg))+len(set(fam))+rank
        for ratios in itertools.product([0., 1.7], repeat=3):
            latent = BasisDraws(dimensions)
            responses = simulate(bg, fam, f, x, beta, scale, ratios, 2*dimensions, latent)
            assert latent.offset == dimensions
            np.testing.assert_allclose(responses.mean(axis=1), x@beta, rtol=0, atol=1e-14)
            centered = responses-(x@beta)[:,None]
            observed = centered@centered.T/(2*dimensions)
            expected = scale*(np.eye(6)+ratios[0]*(bg[:,None]==bg[None,:])+ratios[1]*(fam[:,None]==fam[None,:])+ratios[2]*(f@f.T))
            np.testing.assert_allclose(observed, expected, rtol=1e-13, atol=1e-13)
            maximum = max(maximum, float(np.max(abs(observed-expected))))
            cases += 1
    a = simulate(bg,fam,factor,x,beta,scale,[.3,.4,.5],20,np.random.default_rng(42))
    b = simulate(bg,fam,factor,x,beta,scale,[.3,.4,.5],20,np.random.default_rng(42))
    np.testing.assert_array_equal(a,b)
    invalid = [dict(scale=0),dict(ratios=[-1,0,0]),dict(draws=0),dict(draws=1.5),
               dict(family=np.array(['x','y','x','y','y','y'])),dict(species_factor=np.ones((5,2))),dict(beta=[np.nan,0])]
    base = dict(background=bg,family=fam,species_factor=factor,design=x,beta=beta,scale=scale,ratios=[1,1,1],draws=2)
    for change in invalid:
        try:
            simulate(**(base|change),rng=np.random.default_rng(1))
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid simulation accepted')
    write_json(Path('metadata/matched_working_model_simulator_checks_20260928.json'),dict(
        status='passed_exact_latent_basis_covariance_and_seed_checks',covariance_cases=cases,
        maximum_covariance_error=maximum,invalid_cases_rejected=len(invalid),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/simulate_matched_working_model.py')]},
        scope='Implementation checks only. Basis draws are deterministic quadrature fixtures, not Gaussian samples. Full fitted-data simulations, refitting, interval coverage and model adequacy are not yet evaluated.'))
    print(cases, 'covariance cases;', len(invalid), 'invalid inputs rejected; max error', maximum)


if __name__ == '__main__':
    main()
