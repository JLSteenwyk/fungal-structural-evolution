# Categorical ancestral-state diagnostic

`scripts/ancestral_categorical_diagnostics.py` screens one source-node/extant-
residue anchor across four independently seeded chains. It receives categorical
state codes and constructs a binary indicator for each declared amino-acid,
X or gap state. It never applies a scalar diagnostic directly to arbitrary
integer amino-acid labels. The caller must verify matching model/input,
source-node and anchor identities, independent seeds and burn-in handling.

For nonconstant observed indicators it reuses the existing project's
rank-normalized R-hat, bulk/tail effective sample size and mean Monte Carlo
standard error screen. Its thresholds remain R-hat < 1.01 and both ESS values
at least 400. These diagnostics follow the methods described by
[Vehtari and colleagues](https://arxiv.org/abs/1903.08008); the
[Stan diagnostic reference](https://mc-stan.org/rstan/reference/Rhat.html)
describes rank/folded split diagnostics and approximately 100 effective draws
per chain. Applying them to individual state indicators is this project's
marginal screening design, not a theorem establishing convergence of the
joint ancestral-sequence/alignment posterior.

Constant agreement is reported as requiring review. A state not observed in
any chain has unresolved probability rather than inferred probability zero.
A state that is constant in one chain inherits the conservative review flag
from the existing scalar screen. This is an explicit project diagnostic
policy, not a claim that every constant chain must be scientifically wrong.
We must distinguish strong concentration from failure to explore states.

The report retains each state's four counts and frequencies, all six empirical
pairwise total-variation distances, and every indicator's diagnostic status.
Frequency distances are descriptive: there is no iid categorical test or
confidence interval that ignores chain autocorrelation. Passing observed
indicators is labeled `observed_state_indicator_screens_pass_only`; unseen
states remain explicit and the report never labels a posterior qualified.

The initial saved-sample horizon gives only 75 or 50 retained observations per
chain after the two burn-in cutoffs. The ESS target is not reduced to make
these short traces pass. Adequate mixing may require longer sampling. The
current inference horizon and its pinned scripts are unchanged.

Eight tests passed in the locked ancestral-diagnostic environment: independent
categorical draws with an unseen state; constant agreement; distinct stuck
modes; chain frequency shifts; slow autocorrelation; a rare state absent in
three chains; amino-acid label reordering; and invalid/insufficient draws.
The rare and constant cases cannot return an observed-indicator pass.

```bash
SOFTWARE/ancestral-diagnostics-20260927/bin/python -m unittest discover \
  -s scripts -p test_ancestral_categorical_diagnostics.py
```

This completes the tested diagnostic primitive. Full production integration
still requires provenance-checked four-chain arrays, both cutoffs, preserved
node/anchor denominators, computational planning, and all unresolved states.
No scientific chain quartet was replaced by synthetic replicates, and no
production posterior was assessed or qualified by these fixture results.

## Exact trace reuse

`scripts/ancestral_state_patterns.py` groups only byte-identical complete
chain-by-draw state traces, retaining a pattern identifier for every original
node/anchor coordinate. Chain order, temporal order, state labels and all
coordinates are preserved. Equal frequencies with different draw order remain
different patterns. Identical traces at different anchors do not establish
homology or biological independence: this is computation reuse only.

Five tests passed: full array reconstruction with duplicates; equal frequencies
but different temporal sequences; preserved chain order and state labels;
identity of reused and separately calculated categorical reports at every
fixture coordinate; and rejected malformed inputs. Coordinate multiplicity
must be retained in downstream reporting, even when computation is shared.

The seven completed production chains were each inventoried independently
at both cutoffs. Across those 14 chain/cutoff records, 943,656 node/anchor
records map to 22,129 distinct single-chain temporal patterns. Every input
state reconstructed exactly; no anchor, state, or draw was dropped. The
inventory completed in 1.50 seconds (about 70 MiB measured peak RSS).
Cutoffs overlap, so these totals are not independent biological counts.
Four-chain pattern cardinality can be larger; single-chain compression cannot
be used as a direct estimate of full-grid diagnostic cost. No artificial
quartet was assembled from the seven different input models.

Inventory evidence: `metadata/ancestral_state_pattern_inventory_20260928.json`.
Reproduce with a new receipt path:

```bash
python scripts/inventory_ancestral_state_patterns.py \
  --completion metadata/baliphy_first_anchored_states_completed_20260928.json \
  --output metadata/ancestral_state_pattern_inventory_NEW.json
```

## Verified four-chain array reader

`scripts/read_ancestral_state_quartet.py` binds categorical arrays to the
existing inference quartet verification. It checks distinct chain identifiers
and seeds, one model identity, source configuration/command and artifact
bindings, extraction dispositions and successful attempts, source-audit identity,
biological input alignment, source nodes, tip order and ungapped lengths.
Coordinate metadata must match exactly across chains. The caller must also
verify the producer and extractor plans and their pinned code/input files.

Every array is checked for the expected saved schedule, valid categorical
codes, dimensions and nonnegative unanchored counts. Both burn-in count arrays
are independently recalculated from the full temporal trace. Missing or failed
members return not-ready; no three-chain or mixed-model substitute is made.

Eight handoff tests passed, including the inherited inference provenance tests
and the new full extraction fixture, missing extraction, altered state file,
and invalid counts/schedules/node identities. A complete synthetic quartet
was checked end to end, but these fixtures are not biological posterior samples.

The subsequent readiness check verified the pinned production/extraction plans
and all **405 groups / 1,620 chains**. **Zero groups were ready** in that snapshot.
Evidence: `metadata/ancestral_state_quartet_readiness_20260928.json`.
Production chains continue; the next integration step is the full categorical
report runner, resource accounting and queued application to completed groups.

## Both-cutoff report runner

The report runner now applies the indicator screens to every unique trace at
both cutoffs and retains a pattern-id map for every node/anchor coordinate.
Gzip JSONL records preserve every declared state's per-chain counts,
frequencies and diagnostic status, including unobserved states. Summaries
separately report pattern counts and coordinate-weighted counts, preventing
computational reuse from changing the reported denominator. Neither count
is a count of independent biological sites.

The same report also screens unanchored ancestral residue counts for each
candidate node. This does not assign homology or state identities to those
residues across samples. Full report readback reconstructs every retained
state, verifies all pattern ids/multiplicities and recalculates every saved
per-state count/frequency. The mathematical diagnostic primitive has its
separate tests; report readback is not an independent R-hat implementation.

Three end-to-end tests passed: the provenance-checked synthetic quartet
through both reports with exact coordinate reconstruction; rejection of
duplicate seeds before creating output; and rejection of altered arrays
before creating output. The initial test exposed an unnecessary `psutil`
import from a process helper. The reporter now uses only the already locked
diagnostic dependencies; no environment or running script was changed.

The preparer `scripts/prepare_ancestral_state_quartet_report.py` verifies
producer/extractor plan pins, calls the four-chain reader and writes the bound
array/manifest before invoking `scripts/report_ancestral_state_quartet.py`
under the locked diagnostic Python. Biopython-dependent input checking stays
in the existing analysis environment. Run only after the group has four
verified extracted chains and an appropriate resource plan:

```bash
python scripts/prepare_ancestral_state_quartet_report.py \
  --producer-plan metadata/baliphy_independent_chain_plan_20260927.json \
  --extraction-plan metadata/baliphy_full_anchored_states_plan_20260928_v2.json \
  --group MODEL_INPUT_ID \
  --output results/ancestral/categorical-report-NEW \
  --python SOFTWARE/ancestral-diagnostics-20260927/bin/python
```

The full 405-model categorical queue is now active. At launch, 163 chains had
verified state extractions but no group had all four chains ready. Therefore
no production quartet is yet diagnosed. The controller checks both producer
and extraction process identities and waits for four successful handoffs;
full source verification then precedes each report. Failed or unavailable
groups remain explicit. It never restarts either producer.

The complete input inventory has 24,497,076 node/anchor coordinates per cutoff
(maximum 1,139,968 in one model). Treating every coordinate as a distinct
trajectory at both cutoffs gives 48,994,152 patterns. Applying the largest
measured synthetic per-pattern times gives a planning scenario of 424 CPU-hours
and about 329 GiB of uncompressed report text. These are not hard bounds or
convergence forecasts; quartet compression and data complexity remain unknown.
Single-chain compression was not used to forecast quartet compression.

The queue allows four concurrent reports, four CPU cores, 64 GiB RAM, no swap,
16 GiB address space per process, a 48-hour limit per report, a 512 GiB storage
allowance and at least 768 GiB free disk before launching a report. The active
processing allowance is 12–336 hours, excluding waits for inference chains.
No GPU or paid resources are used.

Eleven queue/handoff tests passed, including an isolated synthetic worker
through both complete reports and verified reuse of its successful attempt,
as well as duplicate-seed rejection and waiting on failed/missing extraction.
Fixtures do not qualify any production posterior.

Reproduce the scheduler with:

```bash
/home/bizon/anaconda3/bin/python scripts/advance_ancestral_categorical_diagnostics.py \
  --plan metadata/ancestral_categorical_queue_plan_20260928.json
```

The immutable plan, launch identity, complete coordinate inventory and tests
are recorded in `metadata/ancestral_categorical_queue_plan_20260928.json`,
`metadata/ancestral_categorical_queue_launch_20260928.json`,
`metadata/ancestral_categorical_full_resources_20260928.json` and
`metadata/ancestral_categorical_queue_checks_20260928.json`.
Outputs go to `results/ancestral/full-categorical-diagnostics-20260928-v1`.
Joint homology/alignment mixing, longer horizons where indicated, root/prior
sensitivity and posterior qualification remain required.
