# BAli-Phy assessment for ancestral uncertainty

The [author guide](https://www.bali-phy.org/README.html) describes joint samples
of ancestral residues and gaps, matched to sampled trees and alignments.
It supports fixed topology or fixed tree analyses. Internal node identifiers
must be interpreted through each matching tree. These capabilities address
uncertainty missing from the current Historian outputs; suitability for our
largest families and MCMC convergence remain untested.

The official [4.3 release](https://github.com/bredelings/BAli-Phy/releases/tag/4.3)
was installed locally using its Ubuntu 24.04 archive. The archive matches the
release SHA-256. Version/help checks pass; no project MCMC was launched.
`install_baliphy_local.py` and its versioned plan reproduce installation without
system changes. Files and help captures are under
`data/software_audits/baliphy-4.3-20260927`.

Before inference, specify priors and indel-rate assumptions, verify all-family
input preservation and capacity, and establish multiple-chain convergence and
sampling diagnostics. Existing Historian checks continue as method sensitivity.

## Initialization checks across the complete selected input grid

`run_baliphy_initialization_grid.py` runs `--test` for all 324 derived-tree
configurations, covering all 78 alignments and retaining all protein copies.
The two floors and all 15 resolutions remain computational sensitivities.
This command initializes a model and exits; it does not sample a posterior.
Resources: two CPU workers, 16 GiB aggregate RAM, no swap, 6 GiB sampled RSS
and 10 minutes per process. The full planning allowance is 0.1–28 hours and
5 GiB output. No GPUs or paid resources are used.

The explicit model arguments are `LG +> ASRV.Gamma(4,alpha=1)`,
`RS07(rate=0.01,meanLength=3)`, constant indel rates, a fixed tree and seed
20260927. The generated Haskell code reveals an additional default: LG
exchangeabilities receive a sampled frequency vector under a symmetric
Dirichlet(1) prior. These frequencies are not the published LG frequencies
or the IQ-TREE empirical estimates. Original internal labels are dropped.
Source data are converted to unaligned character data; source column numbers
must not be used as posterior coordinates. These behaviors are retained and
reported, not silently equated with the preceding model fits.

`audit_baliphy_initializations.py` checks generated model code, input paths,
finite initial scores, prior-plus-likelihood arithmetic, normalization of all
20 frequencies and agreement of total tree length. Its first snapshot verifies
18 jobs and retains 306 pending dispositions. This does not independently
replay likelihoods, verify runtime ancestral identities, establish sample
preservation, measure MCMC mixing or show convergence. Before posterior use,
samples must be linked to the software's internal tree identities and validated
against exact extant sequence identities and candidate descendant sets.

## Same-process sample mapping

The installed Haskell library constructs ancestral names with
`Graph.addAncestralLabel`; its node IDs arise from a process-local ID supply.
A separately loaded tree therefore cannot establish sample identity merely
by matching numerical names. The generated fixed-tree model also omits a tree
logger. For the diagnostic chains, we insert one tree export immediately
before `makeMCMCState`, using `addInternalLabels` and `writeNewick_rooted` on
the actual model tree. No model expressions are changed. The exact original
and modified programs and their hashes are retained for every job.

`run_baliphy_sample_mapping_grid.py` runs 20 iterations for all 324 configurations,
waiting for each verified initialization receipt. Unsuccessful initializations
remain explicit dispositions. This uses two CPU workers, 32 GiB aggregate RAM,
zero swap, 12 GiB sampled RSS and 30 minutes per process, with 20 GiB storage
and a 1–82 hour planning allowance. These are short output-integrity checks,
not final posterior chains, and do not establish MCMC capacity at convergence.

`audit_baliphy_sample_mapping.py` checks the sole tree-export modification,
rooted clades and branch lengths against each specified input tree, all leaf
identities and known residues, and all sample/node associations. The 21 logged
states (iterations 0–20) and three saved alignments (0,10,20) are retained.
Candidate ancestors are identified by exact descendant sets, not numerical
node names. The first four complete chains pass all checks, providing 48
verified candidate-node/sample mappings, with 320 jobs pending in the snapshot.
The degree-two root is retained only conditional on the specified root position.
No early sample is a qualified posterior draw for biological conclusions.

## Effective-input equivalence for future computation reuse

The pinned loader calls `stripGaps` in `mkUnalignedCharacterData`. An audit of
all324 configurations therefore preserves the exact ordered ungapped sequences,
full tree bytes, binary digest and all model/seed options as the computational
input key. It does not sort sequence records or equate different tree files by
approximate distance. FASTA descriptions were checked to equal their identifiers;
only canonical residues,X and gap characters occur. No proteins are removed.

This yields135 distinct input groups:27 groups of four labels and108 groups of
two. The complete324-row mapping retains every alignment, family, tree floor and
resolution label. In85 groups with multiple finished initializations, replacing
only literal input paths gives identical generated program text and exactly
matching initial numerical state lines. Remaining groups still need completed
runtime comparisons; initial agreement is not a mixing or convergence test.

`audit_baliphy_input_equivalence.py` and
`results/ancestral/baliphy-input-equivalence-20260927-v1` preserve the evidence.
Existing diagnostic jobs continue unchanged. Future reuse may map identical
inputs to one computational result per seed, with explicit provenance for every
original label. Different models, priors, trees, seeds and independent chains
must remain separate. A reused result is never an additional independent
replicate. This reduces redundant calculations without reducing family sampling.

## Complete initialization and input-equivalence readback

Both the 324-job initialization producer and its final audit are terminal with
successful exit status. All 324 generated models and finite initial-score
checks pass. The bound artifact hashes were rechecked against the final
completion record. Maximum initialization time was 31.8 seconds and maximum
sampled RSS was 1.05 GiB; these initialization costs do not estimate long-chain
runtime or posterior mixing.

The complete v2 equivalence audit now compares every member of all 135 input
groups. All groups have identical normalized generated programs and initial
numerical-state lines. All 324 labels remain in the configuration mapping;
27 groups contain four labels and 108 contain two. The new script asserts
complete runtime evidence, preserving the earlier partial snapshot separately.
Evidence: `metadata/baliphy_initialization_final_readback_completed_20260927.json`
and `metadata/baliphy_input_equivalence_v2_completed_20260927.json`.
Short-chain runtime/sample integrity checks continue separately. No initialization
result or equivalent-input reuse establishes a qualified posterior ensemble.

A third independent sample-mapping snapshot checks all 36 completed diagnostic
chains, with 288 pending. Every saved alignment, runtime tree, extant sequence
and candidate-node association passes the declared checks, covering 432
candidate-node/sample mappings across 5 families (largest completed
family: 83 proteins). The complete 324-configuration run remains active.
Evidence is recorded in
`metadata/baliphy_sample_mapping_partial_readback_v3_20260927.json`.
These are still 20-step output-integrity diagnostics; the snapshot does not
qualify posterior draws, establish convergence, or cover the unfinished families.


## Measured resource summary for the next sampling stage

`scripts/summarize_baliphy_resources.py` preserves all 324 configurations and
all 135 effective-input groups, joins completed runs to checksum-bound sample
audits, and reports elapsed time, sampled peak RSS and retained artifact sizes.
Failed and pending attempts remain explicit; capped runtimes are censored.
Group tables retain the number of original labels and successful sample audits,
so equivalent inputs cannot be counted as independent chains.

The script was exercised against the frozen 36-chain audit and its total
worker time and pending count were checked separately against the emitted
table. A separate service is queued behind the verified final sample auditor;
it will summarize all 324 dispositions in
`results/ancestral/baliphy-resources-20260927-final-v1`. Its plan and launch
identity are versioned. This is a read-only stage capped at one CPU and 2 GiB,
with no GPU use or paid resources. Measured initialization plus 20-iteration
costs do not estimate iterations required for convergence or production-chain
runtime. Those assessments remain required before qualified posterior use.


## First largest-family short chains independently verified

The v4 sample-mapping audit checks 86 completed configurations and preserves
238 pending dispositions. All 1,032 saved candidate-node/sample mappings pass
the declared integrity checks. This now includes both 622-protein OG0000972
whole-MAFFT configurations, at minimum-edge floors 10⁻⁹ and 10⁻⁷. All known
extant residues, full runtime-tree clades and branches, and candidate identities
were checked across their three saved alignment samples.

The two runs took 1,236.8 and 1,244.7 seconds, with sampled peak RSS of 2.53 and
2.45 GiB, respectively. They therefore fit within the current diagnostic
resource limits; this does not establish the time or memory needed for a
converged chain. Only 20 iterations were run. The corresponding FAMSA and
domain configurations and the remaining full grid continue. Both floor
alternatives share data and are not independent convergence chains.

Evidence: `metadata/baliphy_sample_mapping_partial_readback_v4_20260927.json`
and `metadata/baliphy_resource_partial_v2_completed_20260927.json`. The queued
full audit and all-configuration resource summary remain unchanged.


## Both largest-family alignment labels checked at both floors

The v5 audit passes 88 configurations and all 1,056 saved candidate-node/sample
mappings, retaining 236 pending. All four OG0000972 whole-protein configurations
now pass, including both FAMSA runs. The FAMSA runs completed in 1,230.9 and
1,243.4 seconds. Neither required a memory or time-limit retry.

At each minimum-edge floor, the MAFFT and FAMSA labels have byte-identical
saved alignment samples and byte-identical exported runtime trees. Their
hashes and paths are preserved in
`metadata/baliphy_sample_mapping_partial_readback_v5_20260927.json`. This is
consistent with BAli-Phy's stripping of input gaps and the earlier effective-
input equivalence audit. It is not agreement between independent chains or
support for a particular alignment: the runs share effective data and seed.
All four original dispositions are retained. The complete grid and its final
audits continue; convergence, model/root sensitivity and qualified posterior
samples remain outstanding.

## Installed prior definitions audited before longer sampling

`scripts/audit_baliphy_prior_definitions.py` binds nine installed BAli-Phy 4.3
model/distribution source files by SHA256 and verifies the default expressions.
The Haskell implementation defines LogLaplace as the exponential of a Laplace
variable, with the second parameter its scale. Closed-form quantiles were
checked against SciPy's inverse CDF and by CDF inversion. Results and source
hashes are in `metadata/baliphy_prior_definition_audit_20260927.json`.

| Parameter | Installed default | Median | Central 95% prior interval |
| --- | --- | ---: | ---: |
| Gamma shape alpha | LogLaplace(6,2) | 403.429 | 1.00857–161371.517 |
| RS07 indel rate | LogLaplace(-4,0.707) | 0.0183156 | 0.00220290–0.152283 |

Only 2.489% of the default alpha prior lies below 1. For the unit-mean Gamma
rates, large alpha concentrates rates near 1; this default therefore needs
explicit consideration when assessing among-site rate variation. It must not
be silently substituted for the current diagnostic setting alpha=1. The
diagnostics also fixed indel rate=0.01 and mean indel length=3. Installed RS07
instead defaults to `~ShiftedExponential(10,1)` for meanLength. Installed +F
uses `~SymmetricDirichletOn(letters(@a),1)`; LG supplies exchangeabilities,
while `+> F(LG_freq)` explicitly fixes the published LG frequencies.

Two prospective alpha sensitivity options are recorded: LogLaplace(0,1)
(median 1, central 95% interval 0.05–20) and LogLaplace(0,2) (median 1,
interval 0.0025–400). Neither is selected or launched by this audit. Longer
sampling must explicitly declare all priors, verify the generated model,
rebind the effective-input groups to those settings, and assess independent
chains for mixing and convergence. Source-level defaults and short-chain
capacity checks alone do not qualify ancestral posterior samples.


## Subsequent short-chain snapshot

The v6 independent readback verifies 90 completed configurations and 1,080
saved candidate-node/sample mappings, retaining all 234 pending dispositions.
Every saved sample passes the declared extant-residue, runtime-tree and
candidate-identity checks. Evidence is recorded in
`metadata/baliphy_sample_mapping_partial_readback_v6_20260927.json`. These
remain 20-iteration integrity diagnostics, not converged posterior ensembles.


## Explicit-prior initialization across the full effective-input grid

The next initialization batch covers all 135 distinct ordered-sequence/tree
inputs under each of the three alpha priors above (405 executions). Every
original configuration ID is retained in the representative mapping; aliases
are not counted as independent chains. The preparation script checks input
and mapping hashes against the completed equivalence audit. No inputs are
selected by diagnostic success or biological result.

All three settings explicitly use LG +F, four Gamma categories, the audited
RS07 rate and meanLength priors, fixed input trees, scale 1, and constant
branch indel-rate multipliers. The seed is 20260929. This is `--test`
initialization, not posterior sampling. The package, centered and broad prior
labels describe sensitivity alternatives; no preferred posterior is selected.

The plan reserves two CPU workers, 16 GiB aggregate memory, zero swap and
5 GiB output. Each job has a 600-second and 6-GiB limit; the conservative
runtime envelope is 34 hours including modest overhead, not a convergence
estimate. No GPU or paid infrastructure is used. Plan and launch identity:
`metadata/baliphy_prior_initialization_{plan,launch}_20260927.json`.

The separate auditor checks generated prior expressions and sampled parameter
support, finite initial scores and their arithmetic, normalized frequencies,
fixed input paths, scale and total tree length. Its first frozen snapshot
passes four completed initializations, including all three prior settings for
one input, and retains 401 pending. It does not independently replay the
likelihood or establish convergence. Evidence:
`metadata/baliphy_prior_initialization_partial_readback_20260927.json`.


The complete prior-grid auditor is now queued behind the verified producer
PID and creation time. It requires terminal service success and exit code 0,
then verifies the producer plan and all 405 receipt hashes before auditing.
All failed job dispositions remain explicit; a finished batch does not imply
all model settings passed. The handoff rechecks its pinned files while
waiting and records the final disposition counts only after verifying all
audited artifacts. It is capped at one CPU and 8 GiB, with no GPU use.
Plan and launch: `metadata/baliphy_prior_initialization_final_readback_{plan,launch}_20260927.json`.
The eventual completion record is
`metadata/baliphy_prior_initialization_final_readback_completed_20260927.json`;
its absence means this final audit is not yet complete.


The second independent prior-initialization snapshot checks 140 completed
configurations and retains 265 pending. All completed generated models,
initial score arithmetic, parameter-support and tree-length checks pass.
The per-prior counts and observed initial alpha range are recorded in
`metadata/baliphy_prior_initialization_partial_readback_v2_20260927.json`.
These are initial prior draws, not posterior estimates or independent
convergence chains. The full batch and its final auditor continue.


The v3 prior-initialization readback passes 332 configurations and retains
73 pending. This includes one 622-sequence input under all three explicit
alpha priors, with free RS07 rate/meanLength and frequencies. Each generated
prior expression, initial probability arithmetic, parameter support and
fixed-tree length passes. Per-job runtime and sampled peak RSS for these
largest-input checks are recorded in
`metadata/baliphy_prior_initialization_partial_readback_v3_20260927.json`.
Initialization capacity does not establish longer-chain memory, runtime or
convergence. Other configurations and final readback remain outstanding.


### Completed explicit-prior initialization batch

All 405 configurations (135 effective inputs under each of three alpha
priors) completed with exit code zero. The producer and final auditor both
reached inactive/success with exit code zero. All 405 generated-model,
initial-score arithmetic, parameter-support, frequency-normalization and
tree-length checks passed. Final receipt, source pins and output hashes were
rechecked after service completion. The completion record is
`metadata/baliphy_prior_initialization_final_readback_completed_20260927.json`.

This completes initialization validation across the full effective-input
grid. It does not establish posterior convergence, independent likelihood
validation, or the resource requirements of longer sampling chains. The
separate short-sampling grid and its final resource audit remain pending.


The v7 short-sampling audit verifies all saved samples and candidate-node
identities for 278 of 324 configurations (3,336 candidate-node samples),
retaining 46 pending dispositions. No completed configuration failed these
checks. The v3 resource summary accounts for 19,856.493 worker-seconds
across these 278 runs; the maximum sampled RSS is 2,716,311,552 bytes.
Source and output hashes and resource-table totals were independently
rechecked. These 20-iteration runs include startup overhead and are not
convergence evidence or direct runtime forecasts for longer free-parameter
chains. Full-grid sampling and the queued final audits remain active.
Records: `metadata/baliphy_sample_mapping_partial_readback_v7_20260927.json`
and `metadata/baliphy_resource_partial_v3_completed_20260927.json`.


### Recovery interface review and longer-chain requirements

The installed advanced, expert and developer help and ten installed MCMC
Haskell modules contain no matches for checkpoint, restart, resume,
serialization or restore. This bounded search is reproducible with
`scripts/audit_baliphy_recovery_interface.py`; its receipt is tracked in
`metadata/baliphy_recovery_interface_audit_20260927.json`. It does not
establish that the native implementation or every other interface lacks
checkpoint support. No interrupted-chain continuation has been validated.

The inspected generated program calls `makeMCMCState` before `runMCMC`.
Relaunching that program is therefore not a demonstrated continuation of its
previous chain. The current short-run controller reuses completed,
hash-verified receipts, but an incomplete existing job directory prevents
its `mkdir(exist_ok=False)` path from proceeding. It is not a complete
interrupted-job recovery implementation. Do not modify this pinned controller
while its batch or final audit is active.

The longer-chain controller must use separate immutable attempt directories,
record seed and process identity, and retain every interrupted attempt.
Reuse a successful attempt only after checking its configuration and output
hashes. A rerun starts a new chain with its own burn-in and diagnostics;
never concatenate its samples with an interrupted chain or count reused
seed/input aliases as independent chains. Before replacing any attempt,
verify that its original process and descendants have terminated. Interrupted
attempts remain excluded from posterior qualification unless explicitly
reviewed. Runtime-tree export and iteration-matched alignment readback remain
required for every chain. Native checkpoint continuation would require a
separate interrupted/uninterrupted equivalence check, including random state,
model parameters and saved-sample identities, before use.

These are implementation requirements for the next controller, not a claim
that longer chains, recovery tests or posterior convergence are complete.


### Attempt recovery helper implemented; sampler integration pending

`scripts/ancestral_chain_attempt.py` now implements separate attempt
folders, a nonblocking exclusive lock inherited by the launched process,
configuration/executable/input hash binding, atomic process and completion
records, timeout handling, and process-group checks before reusing output or
starting another attempt. It retains failed and interrupted attempts and
refuses automatic recovery if process identity is missing. Successful exit
is labeled pending scientific validation. Existing completed artifacts are
hash-checked on reuse; altered inputs/configuration or outputs stop reuse.
No samples are concatenated and no native checkpoint is claimed.

Five disposable-process tests pass: reuse and tamper detection, retained
failed attempts, timeout/configuration binding, missing-identity refusal,
and forced parent crash with a live child followed by a separate recovery
attempt. The helper relies on a child retaining its inherited lock descriptor
and remaining in its recorded process group. Actual BAli-Phy behavior must
be qualified before deployment. The surrounding controller must provide
resource limits (including memory), scientific sample audits, stable chain
identifiers and independent seed scheduling. This helper is not yet wired
into a longer-chain batch and does not change the active pinned controller.
Test evidence: `metadata/ancestral_chain_attempt_tests_20260927.json`.


The installed BAli-Phy recovery qualification has now passed using a separate
77-tip OG0002650 input with free parameters and the centered alpha prior.
After its log appeared, the sampler retained inherited lock descriptor 3;
killing its test controller did not permit a duplicate launch. Once that
sampler was terminated, a fresh attempt completed 20 iterations, preserved
all interrupted files, and passed hash-verified reuse without another launch.
All 21 logged states had finite scores and consistent prior/likelihood sums.
The service terminated successfully and all source/output hashes were checked
again independently. Plan, launch and completion records are
`metadata/baliphy_attempt_recovery_{plan,launch,completed}_20260927.json`.

This qualifies that process-recovery path for the installed binary and tested
configuration. It does not resume native MCMC state or qualify ancestral
samples. Memory/resource limits, runtime-tree mapping, independent-chain
seeds and convergence assessment still belong in the longer-chain controller.
No active scientific batch was interrupted by this isolated test.


### Full independent-chain inputs prepared

Prepared 1,620 chain identities: 135 distinct effective inputs × three
explicit alpha priors × four unique seeds. All 324 original configuration
aliases remain attached; aliases never count as separate chains. Each of the
405 generated model programs differs from its audited initialization source
only by an output-relative runtime-tree export immediately before MCMC state
creation. An independent comparison checked every program, input tree and
alignment hash, seed uniqueness and alias membership. No sampling was launched.

Inputs are in `results/ancestral/baliphy-independent-chain-inputs-20260927-v1`;
preparation and check evidence is in
`metadata/baliphy_independent_chain_inputs_completed_20260927.json`.
The complete short-run resource audit, costed iteration horizon, controller
integration, output validation and convergence/extension policy remain required
before launching this grid. Distinct seeds alone do not prove independent
stationary samples or adequate exploration.


### Scalar independent-chain diagnostics implemented

`scripts/ancestral_chain_diagnostics.py` reads four distinct, checksum-bound
logs for one declared model/input identity. It requires unique seed and chain
identifiers and exactly matching iteration schedules; it neither silently
truncates chains nor concatenates attempts. Burn-in is an explicit inclusive
iteration cutoff, applied identically to every chain. Variables must be listed
explicitly in the input manifest. Identity/seed declarations still require
independent provenance checks against the launch and source receipts.

The isolated environment `SOFTWARE/ancestral-diagnostics-20260927` uses
ArviZ 0.22.0; installed versions are frozen in
`environments/ancestral-diagnostics-20260927.lock.txt`. Recreate with a Python
3.10 virtual environment and `pip install -r` that lock file.
The screen reports rank-normalized split/folded R-hat, bulk ESS, tail ESS and
mean Monte Carlo standard error. Initial thresholds are R-hat <1.01 and bulk
and tail ESS >=400 for four chains, following
[Stan diagnostic guidance](https://mc-stan.org/rstan/reference/Rhat.html).
Implementations use the versioned ArviZ
[R-hat](https://python.arviz.org/en/v0.22.0/api/generated/arviz.rhat.html) and
[ESS](https://python.arviz.org/en/v0.22.0/api/generated/arviz.ess.html) APIs.
Any constant chain is explicitly flagged rather than treated as convergence;
truly fixed parameters must be identified from the model and excluded from
monitored stochastic variables deliberately. Nonfinite and very short traces
cannot pass. A passing scalar screen is necessary but not sufficient.

The next controller must include all stochastic model parameters and log
scores, alignment length and indel/substitution summaries. Candidate-node
lengths, positional state probabilities and homology/alignment uncertainty
need separate diagnostics tied to each runtime tree. Compare prespecified
burn-in cutoffs (25% and 50% of each proposed horizon) and report both, rather
than selecting the cutoff that passes. Longer horizons or revised sampling
are required where results disagree; no finite iteration count establishes
convergence. These policies have not yet been applied to a production
independent-chain batch.


### Full short-sampling audit and resource summary complete

All 324 configurations completed successfully and passed the final saved-sample
and candidate-node identity audit: 3,888 candidate-node samples, including all
135 effective-input groups and their aliases. The producer, final auditor and
resource-summary services each terminated inactive/success with exit code zero.
All audit/resource pins and output hashes, all group membership counts, total
worker time and maximum sampled RSS were independently rechecked.

Total measured time across the 324 short runs was 22,901.970 worker-seconds
(6.36 worker-hours); maximum sampled RSS was 2,716,311,552 bytes (2.53 GiB).
These 20-iteration timings include startup and fixed-parameter models; longer
free-parameter chains may differ substantially. Aliases are not independent
chains. Completed integrity checks do not establish convergence or qualify
ancestral structures. Evidence:
`metadata/baliphy_sample_mapping_final_readback_completed_20260927.json` and
`metadata/baliphy_resource_final_completed_20260927.json`.


### Full independent-chain grid launched

The full 1,620-chain grid is now running under
`fungal-baliphy-independent-chains-20260927.service`: 135 effective inputs,
three free-parameter prior settings and four distinct seeds, retaining all
324 original aliases. The first horizon is 1,000 iterations per chain.
It is a sampling checkpoint, not a convergence criterion. The controller
checks every saved alignment and candidate-node identity on successful exit;
capped/failed attempts remain explicit. Scalar and ancestral/alignment mixing
diagnostics and longer horizons remain necessary before sample qualification.

Resource planning uses 10,452.210 seconds summed across effective-input median
20-iteration runtimes. Linear scaling gives 1,742 worker-hours, or 4.54 days
on 16 workers, for this horizon. Startup costs and free-parameter moves make
that an uncertain scheduling baseline, not an ETA or convergence forecast.
The plan records 0.5–4-fold timing sensitivity and per-chain timeouts of at
least one hour (four times the corresponding linear projection where larger).

The batch is limited to 16 CPU equivalents, 192 GiB RAM and zero swap; each
process has a 12-GiB address-space limit and each file a 2-GiB limit. Output
planning allowance is 256 GiB, with at least 3 TiB disk headroom required
before each new attempt. No GPUs or paid resources are used. All 16 initial
workers were verified live with the expected commands, inherited attempt
locks and kernel limits. Initial memory observations do not establish peaks.

The new readback implementation passed fixtures covering all 135 effective
inputs and 1,620 candidate-node samples from the completed short-run batch.
All 1,620 launch configurations, unique seeds and 536 distinct input-file
hashes were independently checked. Plan, launch, prelaunch and startup checks
are tracked under `metadata/baliphy_independent_chain_*_20260927.json`.
The controller and readback code are separate from the completed short-run
scripts; completed short-run evidence is preserved.
