"""Actual frozen working-model refits of dependent Gaussian responses.

This single-model bridge does not admit a production source census or replace
original numerical/timing/fit closures. Conditional normal Wald intervals are
prespecified calibration targets, never asserted to have nominal coverage.
"""
from copy import deepcopy
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.stats import norm

from ancestral_chain_attempt import sha
from full_expanded_model_design_sources import array_digest, digest
from full_weighted_fit_admission import FULL, production_contract
from readback_weighted_shared_entity_candidate import numeric
from reference_measurement_union_sources import verify
from weighted_shared_entity_candidate import READY, candidate, validate_source
from weighted_shared_entity_simulation import GaussianSharedEntitySimulation, coverage_counts, replay_response

SCHEMA = 'weighted-shared-entity-gaussian-refit-v1'
EXPORT_KEYS = {'disposition', 'fit', 'numerical_attempted', 'error_type', 'error_message'}


def original_identity(exported):
    return deepcopy({k: v for k, v in exported.items() if k not in EXPORT_KEYS})


class GaussianRefit:
    def __init__(self, plan_path, expected_plan_sha256, source, rows, identity,
                 expected_identity_sha256, design, original_response, operators,
                 audit, route, diagonal, expected_model_contract):
        self.plan_path = str(plan_path)
        self.plan_sha256 = sha(plan_path)
        if self.plan_sha256 != expected_plan_sha256:
            raise ValueError('Original immutable fit plan hash required')
        self.plan = json.loads(Path(plan_path).read_text())
        verify(self.plan['pins'])
        production_contract(self.plan, dict(logical_cases=FULL['logical_cases'],
            unique_cohorts=FULL['cohorts'], candidate_rows=FULL['candidate_rows'],
            model_setting_rows=FULL['settings']))
        if self.plan['independent_backend'] != 'component_spectral_v1':
            raise ValueError('Original independent backend required')
        self.identity = deepcopy(identity)
        if EXPORT_KEYS & set(identity) or digest(identity) != expected_identity_sha256:
            raise ValueError('Original upstream source identity hash required')
        self.audit, self.route = deepcopy(audit), deepcopy(route)
        d = validate_source(identity, self.audit, self.route, diagonal)
        if identity['source_combined_disposition'] != READY:
            raise ValueError('Unqualified source cannot generate refit evidence')
        rows = np.asarray(rows)
        labels = np.asarray(source['labels'])
        factor = np.asarray(source['factors'][identity['tree']])
        if (rows.ndim != 1 or rows.dtype.kind not in 'iu' or len(rows) != identity['records']
                or len(np.unique(rows)) != len(rows) or np.any(rows < 0) or np.any(rows >= len(labels))
                or labels.ndim != 1 or factor.ndim != 2 or factor.shape[0] != len(labels)):
            raise ValueError('Unique ordered in-bounds original source rows required')
        y = np.asarray(original_response, dtype=float)
        if y.shape != (len(rows),) or not np.isfinite(y).all() or array_digest(y, '<f8') != identity['response_sha256']:
            raise ValueError('Original response does not match source identity')
        self.model = GaussianSharedEntitySimulation(labels[rows], operators, factor[rows], d,
            design, audit['retained_kernel_names'][1:])
        if self.model.input_contract != expected_model_contract:
            raise ValueError('Original upstream design/operator/factor model contract required')
        for name in ['coefficient_columns', 'active_column_indices']:
            if name in identity and len(identity[name]) != self.model.design.shape[1]:
                raise ValueError('Original active coefficient mapping required')
        # Only this selected cohort is copied; the global factor bank is not.
        selected_labels = np.array(labels[rows], copy=True)
        selected_labels.flags.writeable = False
        self.source = dict(labels=selected_labels, factors={identity['tree']: self.model.factor})
        self.rows = np.arange(len(rows))
        self.rows.flags.writeable = False
        self.original_rows_sha256 = array_digest(rows, '<i8')
        self.source_identity_sha256 = expected_identity_sha256
        self.source_contract = digest(dict(identity=identity, ordered_original_rows=self.original_rows_sha256,
            generating_model=self.model.input_contract))
        self.interval = dict(name='conditional_normal_wald_nominal_95_percent',
            nominal_probability=.95, reference='standard_normal',
            critical_value=float(norm.ppf(.975)), covariance='fitted_conditional_beta_covariance',
            variance_estimation_uncertainty_included=False, empirical_coverage_established=False)
        self.environment = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
            platform=platform.platform(), machine=platform.machine())
        self.math_pins = dict(self.plan['pins'])
        for name in ['refit_weighted_shared_entity_simulation.py', 'weighted_shared_entity_simulation.py']:
            p = str(Path(__file__).parent/name)
            self.math_pins[p] = sha(p)
        self.configuration = dict(schema=SCHEMA, original_plan_sha256=self.plan_sha256,
            source_contract=self.source_contract, method=identity['method'],
            optimizer=deepcopy(self.plan['optimizer']), independent_audit=deepcopy(self.plan['independent_audit']),
            independent_backend=self.plan['independent_backend'], interval_definition=self.interval,
            environment=self.environment, math_source_hashes=self.math_pins)
        self.calibration_contract = digest(self.configuration)
        self.norms = np.array([float(z.multiply(z).sum())/self.model.n for z in self.model.incidence.values()]
            + [float(np.sum(self.model.factor*self.model.factor))/self.model.n])/float(np.mean(self.model.diagonal))
        if not np.isfinite(self.norms).all() or np.any(self.norms <= 0):
            raise ValueError('Refitting requires positive original component normalizations')

    def generate(self, beta, scale, ratios, master_seed, scenario_id, replicate):
        self.model.scenario(beta, scale, ratios)
        with np.errstate(over='ignore', invalid='ignore'):
            scaled = np.asarray(ratios)*self.norms
        if not np.isfinite(scaled).all() or np.any(scaled > self.plan['optimizer']['maximum_scaled_variance']):
            raise ValueError('Generating variance lies outside original scaled optimizer box')
        y, generation, _ = self.model.draw(beta, scale, ratios, master_seed,
            self.identity['candidate_id'], scenario_id, replicate)
        return y, generation, self.simulated_identity(generation)

    def simulated_identity(self, generation):
        identity = deepcopy(self.identity)
        identity.update(candidate_id=digest([SCHEMA, self.source_contract, self.calibration_contract, generation]),
            response_sha256=generation['response_sha256'], original_candidate_id=self.identity['candidate_id'],
            original_response_sha256=self.identity['response_sha256'],
            original_source_identity_sha256=self.source_identity_sha256,
            original_source_rows_sha256=self.original_rows_sha256,
            gaussian_generation_sha256=digest(generation), gaussian_model_contract=self.model.input_contract,
            refit_calibration_contract=self.calibration_contract)
        return identity

    def validate_generation(self, response, generation, identity):
        if generation['candidate_id'] != self.identity['candidate_id'] or generation['input_contract'] != self.model.input_contract:
            raise ValueError('Foreign source/model generation')
        y = replay_response(self.model, generation)
        np.testing.assert_array_equal(y, response)
        if identity != self.simulated_identity(generation):
            raise ValueError('Generated source identity does not replay')
        self.model.scenario(**dict(beta=generation['scenario']['beta'], scale=generation['scenario']['scale'],
            ratios=generation['scenario']['variance_ratios']))
        scaled = np.asarray(generation['scenario']['variance_ratios'])*self.norms
        if not np.isfinite(scaled).all() or np.any(scaled > self.plan['optimizer']['maximum_scaled_variance']):
            raise ValueError('Replayed generating variance outside original box')

    def produce(self, response, generation, identity):
        self.validate_generation(response, generation, identity)
        return candidate(self.source, self.plan, self.rows, identity, self.model.design, response,
            self.model.incidence, self.audit, self.route, self.model.diagonal)

    def read(self, response, generation, identity, primary):
        self.validate_generation(response, generation, identity)
        # Export identity corruption is an invariant failure, not a soft numeric review.
        assert original_identity(primary) == identity
        try:
            independent = numeric(self.source, self.plan, self.rows, identity, self.model.design, response,
                primary, self.model.incidence, self.audit, self.route, self.model.diagonal)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            independent = dict(disposition='independent_refit_error_requires_review',
                error_type=type(error).__name__, error_message=str(error), scientific_eligibility=False)
        working = independent['disposition'] == 'independently_audited_working_candidate_pending_inferential_calibration'
        fit = primary['fit']
        status = ('refit_independently_checked_pending_calibration' if working else
            'refit_error_requires_review' if fit is None or 'error_type' in independent else 'refit_requires_review')
        diagnostics = None
        if working:
            beta = np.asarray(fit['beta']); covariance = np.asarray(fit['conditional_beta_covariance'])
            truth = np.asarray(generation['scenario']['beta'])
            assert beta.shape == truth.shape and covariance.shape == (len(beta), len(beta))
            assert fit['coefficient_covariance_is_conditional'] is True
            assert np.isfinite(covariance).all() and np.all(np.diag(covariance) > 0)
            se = np.sqrt(np.diag(covariance)); critical = self.interval['critical_value']
            lower, upper = beta-critical*se, beta+critical*se
            diagnostics = dict(fitted_beta=beta.tolist(), generating_beta=truth.tolist(),
                standard_errors=se.tolist(), standardized_errors=((beta-truth)/se).tolist(),
                nominal_lower=lower.tolist(), nominal_upper=upper.tolist(),
                covers=((lower <= truth)&(truth <= upper)).tolist())
        record = dict(status=status, schema=SCHEMA, candidate_id=generation['candidate_id'],
            scenario_id=generation['scenario_id'], replicate=generation['replicate'],
            input_contract=generation['input_contract'], scenario_contract=generation['scenario_contract'],
            calibration_contract=self.calibration_contract, calibration_configuration=deepcopy(self.configuration),
            generation=deepcopy(generation), simulated_identity=deepcopy(identity), primary=deepcopy(primary),
            independent=independent, nominal_interval_diagnostics=diagnostics,
            fits_computed=1, scientific_eligibility=False, full_production_source_admitted=False,
            nonuniform_weighting_accepted=False, component_variance_attribution_accepted=False,
            scope='Actual frozen single-model working refit; source closure, full generating scenarios/counts, '
                'fixed-sample calibration, adequacy, global dependence/multiple testing and biological acceptance remain separate.')
        if working:record['nominal_interval_covers_generating_beta'] = diagnostics['covers']
        return record

    def replay(self, record):
        if record['calibration_contract'] != self.calibration_contract or record['calibration_configuration'] != self.configuration:
            raise ValueError('Original frozen refit/calibration contract required')
        response = replay_response(self.model, record['generation'])
        fresh = self.read(response, record['generation'], record['simulated_identity'], record['primary'])
        if fresh != record:
            raise ValueError('Original refit/readback/interval/accounting record does not replay')
        return fresh

    def coverage(self, records, replicates):
        # Replay every record before accepting accounting; arbitrary saved status
        # or interval changes must not turn unresolved refits into successes.
        for record in records:self.replay(record)
        result = coverage_counts(records, self.model.design.shape[1], replicates)
        result.update(calibration_contract=self.calibration_contract,
            interval_definition=deepcopy(self.interval), empirical_coverage_accepted=False)
        return result
