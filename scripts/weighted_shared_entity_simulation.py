"""Sparse latent Gaussian draws for exact named shared-entity working models.

Responses preserve D, each named incidence operator and the supplied species
factor. This is a per-model simulation primitive, not inference, a refit, a
global joint-null simulation, or proof that Gaussian assumptions fit the data.
"""
import numpy as np
from scipy import sparse

from full_expanded_model_design_sources import array_digest, digest

SCHEMA='weighted-shared-entity-gaussian-response-v1'


def seed_entropy(master_seed,candidate_id,scenario_id,replicate,input_contract,scenario_contract,component):
    if type(master_seed) is not int or master_seed<0 or type(replicate) is not int or replicate<0:
        raise ValueError('Nonnegative integer master seed and replicate required')
    if any(not isinstance(x,str) or not x for x in [candidate_id,scenario_id,input_contract,scenario_contract,component]):
        raise ValueError('Nonempty source, scenario, contract and component identifiers required')
    payload=digest([SCHEMA,candidate_id,scenario_id,input_contract,scenario_contract,component])
    return [master_seed,replicate,*np.frombuffer(bytes.fromhex(payload),dtype='<u4').tolist()]


class GaussianSharedEntitySimulation:
    """Generate X beta + independent residual/entity/species Gaussian effects.

    With scale s and ratios theta, covariance is
    s*(diag(D)+sum(theta_j*Z_j*Z_j.T)+theta_species*F*F.T).
    Positive D keeps the full covariance nonsingular even at variance boundaries.
    No dense n-by-n covariance or optimizer is used by the generator.
    """
    def __init__(self,labels,incidence,factor,diagonal,design,parameter_names):
        labels=np.array(labels,copy=True)
        if labels.ndim!=1 or not len(labels):raise ValueError('Nonempty block labels required')
        self.n=len(labels);_,blocks=np.unique(labels,return_inverse=True)
        self.diagonal=np.array(diagonal,dtype=float,copy=True)
        if self.diagonal.shape!=(self.n,) or not np.isfinite(self.diagonal).all() or np.any(self.diagonal<=0):
            raise ValueError('Finite positive original residual diagonal required')
        self.factor=np.array(factor,dtype=float,copy=True)
        if self.factor.ndim!=2 or self.factor.shape[0]!=self.n or not np.isfinite(self.factor).all():
            raise ValueError('Finite n-by-r species factor required')
        self.design=np.array(design,dtype=float,copy=True)
        if self.design.ndim!=2 or self.design.shape[0]!=self.n or not 0<self.design.shape[1]<self.n or not np.isfinite(self.design).all():
            raise ValueError('Compatible finite active design required')
        maximum=np.max(abs(self.design),axis=0)
        if np.any(maximum==0):raise ValueError('Inactive zero design column')
        normalized=self.design/maximum;normalized/=np.sqrt(np.sum(normalized*normalized,axis=0))
        singular=np.linalg.svd(normalized,compute_uv=False)
        if singular[-1]<=10*max(self.design.shape)*np.finfo(float).eps*singular[0]:
            raise ValueError('Design rank/boundary requires review')
        self.parameter_names=list(parameter_names)
        if self.parameter_names!=[*incidence,'species'] or len(set(self.parameter_names))!=len(self.parameter_names):
            raise ValueError('Exact unique ordered original component names required')
        self.incidence={};operators={}
        for name,original in incidence.items():
            if not isinstance(name,str) or not name or name in ['species','residual']:
                raise ValueError('Nonempty distinct entity component name required')
            z=sparse.csr_matrix(original,dtype=float,copy=True)
            if z.shape[0]!=self.n or not np.isfinite(z.data).all():raise ValueError('Invalid entity operator')
            z.sum_duplicates();z.eliminate_zeros();z.sort_indices()
            rows,cols=z.nonzero();small=np.full(z.shape[1],self.n,dtype=np.int64);large=np.full(z.shape[1],-1,dtype=np.int64)
            np.minimum.at(small,cols,blocks[rows]);np.maximum.at(large,cols,blocks[rows])
            present=large>=0
            if np.any(small[present]!=large[present]):raise ValueError('Entity spans declared blocks')
            operators[name]=dict(shape=list(z.shape),data=array_digest(z.data,'<f8'),
                indices=array_digest(z.indices,'<i8'),indptr=array_digest(z.indptr,'<i8'))
            for a in [z.data,z.indices,z.indptr]:a.flags.writeable=False
            self.incidence[name]=z
        self.input_description=dict(schema=SCHEMA,records=self.n,block_labels=digest(labels.tolist()),operators=operators,
            species_factor=array_digest(self.factor,'<f8'),residual_diagonal=array_digest(self.diagonal,'<f8'),
            design=array_digest(self.design,'<f8'),parameter_names=self.parameter_names)
        self.input_contract=digest(self.input_description)
        for a in [self.factor,self.diagonal,self.design]:a.flags.writeable=False

    def scenario(self,beta,scale,ratios):
        beta=np.asarray(beta,dtype=float);ratios=np.asarray(ratios,dtype=float)
        if beta.shape!=(self.design.shape[1],) or not np.isfinite(beta).all():raise ValueError('Finite coefficient vector required')
        if ratios.shape!=(len(self.parameter_names),) or not np.isfinite(ratios).all() or np.any(ratios<0):
            raise ValueError('One finite nonnegative variance ratio per original name required')
        if isinstance(scale,(bool,np.bool_)) or not np.isscalar(scale) or not np.isfinite(scale) or scale<=0:raise ValueError('Finite positive residual variance scale required')
        scale=float(scale)
        with np.errstate(over='ignore',invalid='ignore'):
            residual=scale*self.diagonal;components=scale*ratios;mean=self.design@beta
        if not np.isfinite(residual).all() or np.any(residual<=0) or not np.isfinite(components).all() or np.any((components==0)&(ratios>0)) or not np.isfinite(mean).all():
            raise ArithmeticError('Generating moments not numerically representable')
        description=dict(beta=beta.tolist(),scale=scale,variance_ratios=ratios.tolist(),parameter_names=self.parameter_names)
        return mean,residual,components,description

    def draw(self,beta,scale,ratios,master_seed,candidate_id,scenario_id,replicate,keep_latents=False):
        mean,residual,components,description=self.scenario(beta,scale,ratios)
        contract=digest(description);latents={};streams={}
        def normal(name,size):
            entropy=seed_entropy(master_seed,candidate_id,scenario_id,replicate,self.input_contract,contract,name)
            streams[name]=entropy
            latent=np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy))).standard_normal(size)
            if keep_latents:latents[name]=latent
            return latent
        response=mean+np.sqrt(residual)*normal('residual',self.n)
        with np.errstate(over='ignore',invalid='ignore'):
            for j,(name,z) in enumerate(self.incidence.items()):
                response+=np.sqrt(components[j])*(z@normal(name,z.shape[1]))
            response+=np.sqrt(components[-1])*(self.factor@normal('species',self.factor.shape[1]))
        if not np.isfinite(response).all():raise ArithmeticError('Nonfinite generated response requires review')
        record=dict(status='gaussian_response_generated_pending_refit_and_calibration',schema=SCHEMA,
            candidate_id=candidate_id,scenario_id=scenario_id,replicate=replicate,master_seed=master_seed,
            input_contract=self.input_contract,scenario_contract=contract,scenario=description,
            rng='numpy.Generator.PCG64.SeedSequence',component_seed_entropy=streams,
            response_sha256=array_digest(response,'<f8'),scientific_eligibility=False,fits_computed=0,
            scope='Per-model Gaussian working response, not calibrated inference or global joint-null/multiple-testing simulation. Original source admission, generating scenarios, full refits, review accounting, adequacy and lineage/selection uncertainty remain required.')
        return response,record,latents


def coverage_counts(records,coefficients,replicates):
    """Retain every prespecified replicate, including pending/failed refits.

    Fixed-sample exact Monte Carlo envelopes use the previously qualified
    binomial implementation. They are not coefficient intervals or optional-
    stopping guarantees and cannot combine different generating models.
    """
    from matched_calibration_intervals import exact_coverage_bounds
    if type(coefficients) is not int or coefficients<1 or type(replicates) is not int or replicates<1:
        raise ValueError('Positive coefficient and prespecified replicate counts required')
    if len(records)!=replicates or any(type(r['replicate']) is not int for r in records) or sorted(r['replicate'] for r in records)!=list(range(replicates)):
        raise ValueError('Exactly every prespecified replicate required')
    keys={(r['candidate_id'],r['scenario_id'],r['input_contract'],r['scenario_contract']) for r in records}
    if len(keys)!=1:raise ValueError('One exact generating model/scenario per accounting group required')
    covered=np.zeros(coefficients,dtype=int);qualified=0
    for record in records:
        if record['scientific_eligibility'] is not False:raise ValueError('No simulation eligibility promotion')
        if record['status']=='refit_independently_checked_pending_calibration':
            values=record['nominal_interval_covers_generating_beta']
            if len(values)!=coefficients or any(type(v) is not bool for v in values):raise ValueError('Exact boolean coverage vector required')
            covered+=values;qualified+=1
        elif record['status'] not in ['gaussian_response_generated_pending_refit_and_calibration','refit_requires_review','refit_error_requires_review']:
            raise ValueError('Unknown refit state requires explicit review accounting')
    unresolved=replicates-qualified
    return dict(attempted=replicates,independently_checked=qualified,unresolved=unresolved,
        covered_qualified=covered.tolist(),marginal_monte_carlo_intervals=[exact_coverage_bounds(int(k),unresolved,replicates) for k in covered],
        scientific_eligibility=False,scope='Fixed prespecified independent replicates of one generating model; unresolved outcomes remain in denominator and exact interval envelope. Not biological coefficient intervals, simultaneous coverage, or calibration acceptance.')


def replay_response(model,record):
    """Regenerate the entire response/stream contract, without an optimizer."""
    scenario=record['scenario']
    response,fresh,_=model.draw(scenario['beta'],scenario['scale'],scenario['variance_ratios'],
        record['master_seed'],record['candidate_id'],record['scenario_id'],record['replicate'])
    if record!=fresh:raise ValueError('Original response/stream/scenario record does not replay exactly')
    return response
