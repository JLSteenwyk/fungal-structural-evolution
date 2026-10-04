# Dependent Gaussian responses for future full-model calibration

`weighted_shared_entity_simulation.py` now generates responses from the exact
named shared-entity working covariance used by the weighted fitting code.
It preserves positive residual diagonal D, sparse gene/family/pair operators,
signed/coalesced loadings, the supplied phylogenetic factor F and original X.
It uses independent latent Gaussian draws with mean X beta and covariance
`scale * (diag(D) + sum(ratio_j * Z_j Z_j.T) + ratio_species * F F.T)`.
Scale is residual **variance**, not standard deviation. Zero variance
components remain explicit; nonuniform D and target identity stay distinct.
Generation avoids dense observation-by-observation covariance matrices.

This is a **per-model simulation primitive**, not a completed calibration,
native refit, biological pilot or global joint-null simulation. Full original
source/fit admission, generating scenarios, replicate counts, complete refits,
interval/test definitions, selection/phylogenetic uncertainty, shared effects
across different models and multiple-testing control remain to implement or
qualify. Independently generating marginal models does not preserve every
cross-model dependency needed for a simultaneous null experiment.

## Reproducible streams and retained uncertainty

Streams depend on master seed, original candidate identity, scenario identity,
replicate index, input/scenario hashes and named component. Worker order and
scheduling do not select the stream. Every record retains these inputs,
component entropy, original generating parameters and response hash; replay
regenerates the entire record/response and rejects changed fields. Model inputs
are copied and validated; block-spanning entities, invalid/underflowed variances,
incompatible names and rank-boundary designs require review.

Coverage accounting requires exactly every prespecified replicate of one
candidate/input/scenario. Pending, failed and reviewed refits remain in the
denominator. It reuses the previously qualified fixed-sample binomial interval
envelope over every unresolved completion. These intervals describe Monte Carlo
uncertainty, not biological coefficient uncertainty or optional-stopping validity.
Independent refit statuses in these accounting tests are **synthetic fixture
states**; no optimizer or biological refit ran.

The exact-interval reference is
[SciPy 1.15.3 binomtest](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.binomtest.html).
The working Gaussian distribution is defined by its mean/covariance as in
[SciPy 1.15.3 multivariate_normal](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.multivariate_normal.html).
These software references do not establish model adequacy for fungal data.

## Complete declared software checks

Checks reconstruct all 64 previously saved synthetic weighted models,
covering every common/pair-exception, signed/unsigned, four-control,
two-outcome and ML/REML combination, with q4/q5/q6 bases. They retain the
original full production header of 20,832,000 potential candidates and
49,766,400 setting links; software fixtures cannot reduce or admit that scope.
The original model, numerical providers and qualified parent bytes are unchanged.

- 192 dense covariance comparisons across zero, positive and sparse-boundary
  component scenarios; maximum absolute discrepancy 5.33e-15.
- 576 actual generated responses compared with independent explicit latent
  column reconstruction and exact serialized replay; maximum discrepancy
  3.55e-15. Reversed scheduling preserves every response.
- 16 seeded Gaussian moment checks, 512 draws each, supplement the exact
  algebra. These 8,192 draws are software smoke checks, not coefficient coverage
  estimates or production simulation.
- Signed duplicate cancellation, explicit zero incidence and zero-rank species
  factors remain valid; shared signed effects produce negative covariances.
- 29 malformed input/scenario cases, 11 altered replicate accounting cases and
  seven altered response/seed/stream records are rejected.

The original bounded software wait **62497 exited zero**, with exact wrapper
PID 89685/create/command and invocation `6162a20feabe48ad895c0427c2c2a2b3`,
two original wrapper messages and original manager start/completion verified.
There are 500 transport bindings. Native software CPU was 9.72 seconds and
child peak RSS 143,106,048 bytes. The unit used two CPUs, 16 GiB, zero swap,
one BLAS thread, 12-GiB native address-space, 3,600 CPU-second and 256-MiB
file limits; output allowance 2 GiB, wall budget 7,200 seconds. These are
software costs and limits, not full-data refit estimates.

- [Validation](../metadata/weighted_shared_entity_simulation_software_validation_20261004_v1.json).
- [Original execution](../metadata/weighted_shared_entity_simulation_software_execution_20261004_v1.json).
- [Actual wait and original journals](../metadata/weighted_shared_entity_simulation_software_transport_20261004_v1.json).
- [Entire original terminal journal payload](../metadata/weighted_shared_entity_simulation_original_terminal_payload_20261004_v1.json).
- [Recorded installed environment](../metadata/weighted_shared_entity_simulation_environment_20261004_v1.json).
- [Data dictionary](../metadata/weighted_shared_entity_simulation_data_dictionary_20261004_v1.tsv).

The installed environment was observed after execution: Python 3.10,
NumPy 2.2.6 and SciPy 1.15.3, including native RNG/array/statistics artifact
hashes. The original checker did not serialize package versions per response.
Exact replay is tested in this environment; arbitrary version/platform
portability is not asserted. NumPy documents limitations on cross-environment
random-stream compatibility in its
[compatibility policy](https://numpy.org/doc/2.2/reference/random/compatibility.html).

Software artifacts remain outside Git at
`data/software_audits/weighted-shared-entity-simulation-20261004-v1/`.
The checker regenerates fixtures under fresh output/receipt namespaces;
qualified roots and original receipts must not be overwritten. Reproduce under
the original bounded wrapper/resource environment with a fresh namespace:

```bash
python scripts/check_weighted_shared_entity_simulation.py \
  --output data/software_audits/weighted-simulation-NEW \
  --receipt metadata/weighted-simulation-NEW.json
```

No production simulation, native refit, calibrated effect, accepted weighting
or component attribution is complete. The full atlas, accepted phylogenetic/
reconciliation/dating framework and all eight biological aims remain required.
