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
