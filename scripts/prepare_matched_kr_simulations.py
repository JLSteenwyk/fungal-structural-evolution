"""Freeze synthetic method-validation designs, fixed replicate counts and budget."""
import itertools
from pathlib import Path
import numpy as np
from ancestral_chain_attempt import sha,write_json


def main():
    root=Path('results/model_validation/matched-kr-simulation-inputs-20260928-v1')
    root.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(20260928);cases=[];pins={}
    for n,groups,per_family,rank,p in [(48,16,4,4,2),(240,60,5,16,5)]:
        bg=np.repeat(np.arange(groups),n//groups);fam=bg//per_family
        arrays=dict(background=bg,family=fam,factor=rng.normal(size=(n,rank))/np.sqrt(rank),
                    design=np.column_stack([np.ones(n),rng.normal(size=(n,p-1))]),
                    beta=np.array([.1,-.4,.2,-.1,.3])[:p])
        path=root/f'design-{n}.npz';np.savez_compressed(path,**arrays)
        with np.load(path,allow_pickle=False) as saved:
            for key,value in arrays.items():np.testing.assert_array_equal(saved[key],value)
        pins[str(path)]=sha(path)
        for index,ratios in enumerate(itertools.product([0.,.7],repeat=3)):
            cases.append(dict(id=f'n{n}-variance-{index}',design=str(path),records=n,
                              ratios=list(ratios),scale=1.3,coefficients=p))
    names=['prepare_matched_kr_simulations','run_matched_kr_simulations','calibrate_matched_kr',
           'refit_matched_simulation','simulate_matched_working_model','refine_matched_reml_analytic',
           'matched_reml_gradient','cached_matched_likelihood','matched_mixed_covariance',
           'matched_covariance_information','matched_covariance_information_fast','matched_kr_covariance',
           'matched_kr_scalar','matched_calibration_intervals','ancestral_chain_attempt']
    for name in names:
        path=Path('scripts')/(name+'.py');pins[str(path)]=sha(path)
    plan=dict(output='results/model_validation/matched-kr-simulations-20260928-v1',
        cases=cases,replicates=999,master_seed=20260928,pins=pins,
        resources=dict(workers=4,cpus=4,memory_gib=8,swap_gib=0,storage_gib=4,
                       minimum_free_gib=100,planning_hours=[1,24],gpu=False,paid_cost=0),
        resource_basis='15984 refits with24 optimization candidates each; approximate 0.5–20seconds '
                       'per attempt gives0.6–22.2hours at4workers before overhead. Planning range1–24hours; '
                       'service limit24hours.4GiB allowance versus about0.8GiB at50KiB per disposition.',
        statistical_scope='Two frozen Gaussian designs times eight zero/moderate variance configurations. '
                          '999 independent seeded responses per generating model, no optional stopping. '
                          'At true coverage0.95, binomial MCSE is0.006896. Marginal interval envelopes '
                          'retain unresolved attempts; comparisons across cases are descriptive. '
                          'This method-validation experiment neither replaces the full fungal cohort '
                          'nor establishes coverage for every real design or misspecified model.')
    target=Path('metadata/matched_kr_simulations_plan_20260928.json');assert not target.exists()
    write_json(target,plan);print('Frozen',len(cases)*plan['replicates'],'attempts; plan',sha(target))


if __name__=='__main__':main()
