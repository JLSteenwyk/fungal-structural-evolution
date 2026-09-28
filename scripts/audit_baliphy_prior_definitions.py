#!/usr/bin/env python3
"""Bind installed model priors and calculate their implied quantiles."""
import hashlib,json,math
from pathlib import Path
from scipy.stats import laplace


def main():
    root=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/lib/bali-phy')
    files=['bindings/models/RS07.json','bindings/models/ASRV/Gamma.json','bindings/models/F.json','bindings/models/LG.json','bindings/distributions/LogLaplace.json','bindings/distributions/ShiftedExponential.json','haskell/Probability/Distribution/Transform.hs','haskell/Probability/Distribution/Laplace.hs','haskell/Probability/Distribution/Exponential.hs']
    pins={str(root/name):hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files}
    models={name:json.loads((root/'bindings/models'/name).read_text()) for name in ['RS07.json','ASRV/Gamma.json','F.json']}
    defaults={name:{a['name']:a['default_value'] for a in model['args'] if 'default_value' in a} for name,model in models.items()}
    assert defaults['RS07.json']['rate']=='~LogLaplace(-4,0.707)'
    assert defaults['RS07.json']['meanLength']=='~ShiftedExponential(10,1)'
    assert defaults['ASRV/Gamma.json']['alpha']=='~LogLaplace(6,2)'
    assert defaults['F.json']['pi']=='~SymmetricDirichletOn(letters(@a),1)'
    transform=(root/'haskell/Probability/Distribution/Transform.hs').read_text()
    assert 'logLaplace m s = ExpTransform $ laplace m s' in transform
    assert 'sample (ExpTransform dist) = exp <$> sample dist' in transform
    assert '2*s^2' in (root/'haskell/Probability/Distribution/Laplace.hs').read_text()
    def quantiles(location,scale):
        values={}
        for probability in [.025,.5,.975]:
            logq=location+(scale*math.log(2*probability) if probability<=.5 else -scale*math.log(2*(1-probability)))
            value=math.exp(logq)
            assert math.isclose(value,math.exp(laplace.ppf(probability,loc=location,scale=scale)),rel_tol=1e-12)
            assert math.isclose(laplace.cdf(math.log(value),loc=location,scale=scale),probability,abs_tol=1e-12)
            values[str(probability)]=value
        return values
    result=dict(status='complete_installed_prior_semantics_and_quantiles',source_hashes=pins,installed_defaults=defaults,
        default_gamma_shape=dict(distribution='LogLaplace(6,2)',quantiles=quantiles(6,2),probability_below_one=float(laplace.cdf(0,loc=6,scale=2))),
        default_indel_rate=dict(distribution='LogLaplace(-4,0.707)',quantiles=quantiles(-4,.707)),
        current_diagnostic_fixed_values=dict(gamma_shape=1,indel_rate=.01,mean_indel_length=3),
        candidate_gamma_prior_sensitivities=[dict(distribution='LogLaplace(0,1)',quantiles=quantiles(0,1)),dict(distribution='LogLaplace(0,2)',quantiles=quantiles(0,2))],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Installed 4.3 source audit and mathematical prior summaries, not generated-model verification or fitted posterior evidence. Candidate priors are prospective sensitivity options; no production chain or prior selection is implied. Changing models/priors requires new initialization, input binding and convergence assessment.')
    Path('metadata/baliphy_prior_definition_audit_20260927.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['default_gamma_shape','default_indel_rate','candidate_gamma_prior_sensitivities']}))


if __name__=='__main__':main()
