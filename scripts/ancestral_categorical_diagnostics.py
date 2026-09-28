#!/usr/bin/env python3
"""Four-chain state-indicator screens; observed agreement is not convergence."""
import itertools
import numpy as np
from ancestral_chain_diagnostics import diagnose


def diagnose_states(values, alphabet):
    """Screen one node/anchor without imposing an order on state labels.

Input axes are chain, retained saved iteration. Chain provenance, identical
node/anchor identities, seed independence and burn-in must be checked by the
caller. Every declared state is reported, including states never observed.
"""
    values=np.asarray(values)
    if values.ndim!=2 or values.shape[0]!=4 or values.shape[1]==0:
        raise ValueError('Exactly four nonempty equally sized chains required')
    if not alphabet or len(set(alphabet))!=len(alphabet):
        raise ValueError('Distinct declared state labels required')
    if not np.issubdtype(values.dtype,np.integer) or np.any(values<0) or np.any(values>=len(alphabet)):
        raise ValueError('State codes must be valid integer alphabet indices')
    counts=np.array([[np.count_nonzero(chain==s) for s in range(len(alphabet))] for chain in values])
    frequencies=counts/values.shape[1]
    comparisons=[dict(chains=[a,b],total_variation=float(np.abs(frequencies[a]-frequencies[b]).sum()/2))
                 for a,b in itertools.combinations(range(4),2)]
    indicators={}
    for state,label in enumerate(alphabet):
        c=counts[:,state]
        if c.sum()==0:
            screen=dict(status='state_not_observed_probability_unresolved',draws_per_chain=values.shape[1])
        elif np.all(c==values.shape[1]):
            screen=dict(status='state_constant_in_all_chains_requires_review',draws_per_chain=values.shape[1])
        else:
            screen=diagnose((values==state).astype(float))
        indicators[label]=dict(counts_per_chain=c.tolist(),frequency_per_chain=frequencies[:,state].tolist(),screen=screen)
    observed=[r for r in indicators.values() if sum(r['counts_per_chain'])]
    if len(observed)==1:
        status='no_observed_state_variation_requires_review'
    elif all(r['screen']['status']=='passes_scalar_screen_only' for r in observed):
        status='observed_state_indicator_screens_pass_only'
    else:
        status='categorical_mixing_requires_review'
    return dict(status=status,draws_per_chain=values.shape[1],observed_states=len(observed),
        unobserved_states=[s for s,r in indicators.items() if not sum(r['counts_per_chain'])],
        indicators=indicators,pairwise_empirical_total_variation=comparisons,
        scope='Indicators of declared categorical states; no numeric ordering of amino-acid labels. '
        'Empirical frequency distances are descriptive, without iid tests or uncertainty claims. '
        'Unseen states remain unresolved. No joint alignment/sequence posterior qualification.')
