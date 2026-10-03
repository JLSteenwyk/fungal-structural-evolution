# Full-grid model-preserving ancestral initialization

The new [startup preflight](../metadata/baliphy_reference_preflight_plan_20261003.json)
initializes all **1,620 roles: 135 effective inputs × three original priors ×
four fresh seeds**, retaining all 324 configuration aliases. These are the
selected ancestral case families across their complete original protein/taxon
sampling, rather than the entire structural atlas. Native execution uses
`--test`; it does not generate a new MCMC horizon. The two failed earlier
chains, inadequate posterior mixing, root/model uncertainty and all eight
biological aims remain unresolved.

## October 3: complete native startup and corrected reader closure

All 1,620 original startup-only attempts and their first serialized reader
have completed [full provenance closure](../metadata/baliphy_reference_preflight_completed_20261003.json).
The first parser accepted 1,548 and retained 72 invalid-output dispositions.
The [separate full footer-aware replay](../metadata/baliphy_reference_startup_footer_completed_20261003.json)
has completed with all 1,620 startups valid, 405 complete quartets and zero
failed startup dispositions. All 72 differences are recognized native timing
footers; no native attempt or prior source/output was changed or rerun.
Corrected closure binds 13,960 sources/artifacts and both original journals.

The [fresh compact completion observation](../metadata/ancestral_qualification_completion_execution_checkpoint_20261003.json)
checks all 982 footer-stage pins, archived summaries/receipt bindings and all
three original terminal handles. Earlier dated partial observations below
retain their original scope. The full short sampler has now started after
this gate; its [workflow](baliphy-reference-sampler-qualification-20261003.md)
and [resource observation](baliphy-sampler-resource-observation-20261003.md)
are separate stages. Startup fidelity does not qualify their outputs, prove
long-run memory safety, resolve historical allocation failures or establish
posterior mixing. The longer posterior proposal remains unlaunched.

## What changes and what is preserved

The installed [PhyloAlignment module at the pinned upstream commit](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/haskell/Probability/Distribution/PhyloAlignment.hs)
initializes branch alignments from left-aligned tip sequences and median
internal sequence lengths. It stores fixed lengths only for labeled tips;
ancestral lengths are derived from the incident pairwise alignments. Large
early alignment expansion in the two earlier failed attempts motivated testing
an alternative initial state. This does not establish the cause of either
failure or prove that initialization resolves it.

[The new initializer](../scripts/baliphy_reference_initialization.py) obtains
the original node count and fixed-tip map from `sample_alignment`, then
substitutes the branch alignments constructed from the supplied MAFFT/FAMSA
matrix. It retains the original distribution object in `RanDistribution3`,
`alignment_effect`, `triggeredModifiableAlignment` and property lookup.
The [native Random interpreter](https://github.com/bredelings/BAli-Phy/blob/80b0402eed0157f31ecb57e0efc34c03ed83050c/haskell/Probability/Random.hs)
therefore registers the original density and transition effects. Both original
alignment moves remain enabled. Internal lengths are not fixed to their
starting values.

The matrix is used only to choose the initial state. Observed data remain the
same unaligned amino-acid sequences in the original `phyloCTMC` call. The LG
model, sampled amino-acid frequencies, four gamma categories, three gamma-shape
priors, RS07 indel priors, fixed tree, branch parameters, scale and logging
intervals remain unchanged. Each of six explicit source substitutions must
occur once; reversing them must recover the complete original program byte
for byte. No original model source, native installation or running plan was
edited. All 405 transformed programs are separate outputs.

The [version-matched source audit](../metadata/baliphy_reference_initialization_source_audit_20261003.json)
checks installed PhyloAlignment, Random, Bio.Alignment, Bio.Sequence and MCMC
module bytes against upstream commit `80b0402eed0157f31ecb57e0efc34c03ed83050c`.
The executable is separately hashed and exercised; this is not an independent
rebuild of its binary.

## Checks and full-scope resource plan

The [native software fixtures](../metadata/baliphy_reference_initialization_software_validation_20261003_v3.json)
run 12 revised and 12 original startups on a synthetic four-tip gapped matrix.
All matching-seed stochastic parameter draws agree exactly. Revised startups
preserve extant residue homology, four fixed tip lengths and three free
ancestral lengths. Three 20-iteration software fixtures, one per prior, confirm
that alignment transitions actually move while retaining observed residues.
Test-only audit fields are absent from ordinary inference logs. These fixtures
are software checks, not a fungal pilot or biological posterior analysis.

The annotated alignment log density is compared with a degree-aware evaluation:

`log(branch probability product) + Σ_internal (1 − degree(v)) log(length prior(L_v))`.

The alternative legacy `alignment_pr` helper uses a squared internal correction;
it differs at a degree-two root. An early prototype exposed that distinction.
The verification was changed to the formula used by the actual annotated
distribution; no production density was changed. The native formula comparison
uses absolute tolerance 1e-8 and relative tolerance 1e-12. Score arithmetic uses
absolute tolerance 1e-7 and relative tolerance 1e-12. Existing scientific
convergence thresholds and earlier diagnostic comparison tolerances remain
unchanged.

The native checks reject 23 malformed matrices, missing transformation anchors
and altered audit exports. [The complete generator check](../metadata/baliphy_reference_preflight_grid_validation_20261003.json)
constructs every real program and all 1,620 configurations, rejects ten altered
designs plus a missing disposition, and retains two explicitly simulated failed
roles and unresolved quartets. Simulation tests accounting only. The actual
native attempt decoder additionally [passes all 12 saved native fixtures](../metadata/baliphy_reference_native_attempt_decoder_validation_20261003.json).

Two earlier software-fixture runs failed on the report's `prior` key and the
MCON header expectation. Both were corrected before the final stable version3
validation. Their output roots remain preserved. The prototype ambiguity in
property imports was corrected with qualified `Effect.getProperties`.

| Enforced control | Startup execution |
| --- | --- |
| Total CPU quota / concurrent workers | 2 / 2 |
| Whole-group memory / swap | 24 GiB / zero |
| Per-process address space | 12 GiB |
| Per-process CPU / wall timeout | 600 / 900 seconds |
| Per-file size cap | 256 MiB |
| BLAS/OpenMP threads | 1 |
| Planned output allowance / minimum free disk | 16 / 64 GiB |
| GPU / new charges | None / zero |

Historical 405-startup elapsed worker duration was 4,875.26 seconds. Four-role
linear scaling gives 19,501.03 seconds, with a planning sensitivity range of
9,750.51–78,004.10 worker seconds. The different initial state and current
contention make these uncalibrated scenarios. The summed wall timeout ceilings
are 405 worker hours, not a forecast or convergence estimate. Output allowances
are planning estimates; per-file, CPU, address-space and cgroup caps are
enforced. Available memory and disk were inspected before launch.

## Timing footer integration and complete reread

The native CLI appends a five-line `Work:` timing summary after sufficiently
long startups. The first full-grid decoder required the entire stdout to be
JSON and retains those cases as invalid reader dispositions. Native exit-zero
outputs, sources, seeds, models and original classifications remain unchanged.
They are not silently relabeled inside the original running stage.

A [separate reader version](../scripts/baliphy_reference_startup_readback_v2.py)
accepts one initial JSON model followed by whitespace or this exact footer:
start/end dates and nonnegative reported elapsed/CPU seconds. Unknown text,
second JSON records and incomplete or malformed timing summaries fail. It
rechecks the same model source, configuration, actual native command/process,
receipt hashes, runtime tree, residue homology, tip/free-ancestor representation
and degree-aware density. It never launches another native job.

[Actual-output regression](../metadata/baliphy_reference_startup_footer_validation_20261003.json)
passed for 24 complete-grid startups, including 12 native timing footers, and
rejected seven malformed stdout cases. [Complete serialization contracts](../metadata/baliphy_reference_footer_workflow_validation_20261003.json)
exercise all 1,620 metadata identities with explicitly artificial status rows,
retaining two failures and rejecting eight altered exports and two completed
restarts. This software workflow mocks native reconstruction and does not
substitute for the full real-data readback.

The [complete follow-up](../metadata/baliphy_reference_startup_footer_plan_20261003.json)
is queued after the original native producer, serialized reader and provenance
closure. Its producer and reader will each reread all 1,620 native outcomes.
Closure requires every source/artifact binding and both original completion
journals. Original classifications remain linked to any new classifications.
All three follow-up services retain two-CPU/24-GiB/no-swap/BLAS-one caps. No
native attempt is repeated for this integration repair.

At 05:21 UTC (01:21 EDT), [the actual snapshot replay](../metadata/baliphy_reference_startup_footer_execution_checkpoint_20261003.json)
checked every available 136 startups: all passed, including 24 timing footers
previously rejected by the first decoder. Full native startup/readback/closure
is still pending. At 05:22 UTC [the original native observer](../metadata/baliphy_reference_preflight_execution_checkpoint_20261003_v3.json)
verified 162 unclosed checkpoints, all 974 frozen pins, exact original
controllers 2786767/2786771/2786776, live cgroups and native worker CPU/address-
space/file limits. Follow-up controllers 2792240/2792244/2792248 were separately
verified as the original queued handles.

## Numerical interpretation and remaining decisions

The earlier full census retains 336 initialization observations of infinite
gamma shape, all before burn-in. The pinned native implementation explicitly
treats positive infinite gamma shape as its equal-rate limit. These values are
supported limit behavior requiring interpretation, not by themselves evidence
of invalid retained samples. The measured counts and frozen census remain
unchanged. See [the earlier source review](baliphy-full-first-horizon-diagnostics-20261002.md).

Startup fidelity and better initial geometry cannot qualify a posterior or
establish safe memory for subsequent alignment proposals. Longer fresh chains
still require a full resource/output plan, retained-failure accounting and all
scalar, length, alignment and categorical mixing screens. Root, predictor,
model-adequacy and circularity controls remain separate. The proposed
10,000-iteration posterior grid remains unlaunched; GPU prediction stays
paused. All eight biological aims remain incomplete.


At 05:27 UTC, [a later complete available-startup snapshot](../metadata/baliphy_reference_startup_footer_execution_checkpoint_20261003_v2.json)
checked 234 native outcomes: all passed, including 24 timing footers.
The earlier [full native residue replay](independent-native-alignment-replay-20261002.md#october-3-complete-serialized-and-provenance-closure)
also completed its full 53,331-binding/two-journal provenance closure.
Startup qualification remains in progress; neither integrity stage resolves
posterior mixing or either historical allocation failure.


## Reproducibility and inspecting existing execution

The [CPU environment](../environments/baliphy-reference-startup-20261003.yml)
records the Python dependencies. Native software and every model/input hash
are pinned in the execution plans. Generated models and large native artifacts
remain outside Git; their paths, hashes and transformation records are
bound by the plans and completed archives. The repository retains the
initializer, native fixture checks, full generator, attempt runner, both
readers, launchers and exact-handle observers.

Inspect the already running jobs without launching or restarting inference:

```bash
python scripts/record_baliphy_reference_preflight_checkpoint.py \
  --output NEW_REFERENCE_STARTUP_CHECKPOINT.json
python scripts/record_baliphy_reference_footer_checkpoint.py \
  --output NEW_REFERENCE_FOOTER_CHECKPOINT.json
```

Output names must be new. Complete original native startup closure and the
queued complete footer-aware producer/readback/provenance closure are required
before treating the whole grid as startup-qualified. All simulation and native
software-fixture evidence remains separate from that production claim.

[The development audit](../metadata/baliphy_reference_initialization_development_audit_20261003.json)
preserves prototype/fixture output hashes and the corrected report/header/
timing-grammar expectations separately from the stable software gates.

At 05:32 UTC, [the latest complete available-startup snapshot](../metadata/baliphy_reference_startup_footer_execution_checkpoint_20261003_v3.json)
checked all 312 available records: 312 passed, including 24 timing-footer
reclassifications. The [same-original native observation](../metadata/baliphy_reference_preflight_execution_checkpoint_20261003_v4.json)
rechecked both frozen inputs and worker limits. The full 1,620-role startup
qualification remains incomplete.
