# Full ancestral sampling resource and numerical survey

The [complete survey](../metadata/baliphy_horizon_resources_completed_20261003.json)
accounts for all 1,620 original chain identities, 405 model quartets and 135
effective inputs, including all 324 original configuration aliases. It measures
all 1,620 initial attempts plus the three selected recovery attempts. The
selected set retains 1,618 intact chains and both failures; all six intact
chains in the two unresolved quartets remain included.

The producer, full serialized reconstruction and provenance closure completed
with 8,690 source/artifact bindings and both original producer/reader journals.
The [final execution observation](../metadata/baliphy_horizon_resources_execution_checkpoint_20261003_v2.json)
verified successful termination of controllers 2757537, 2757541 and 2757546,
the unchanged 25-pin plan and the completion archive hash. No native sampling,
GPU prediction or paid resource was launched.

## Measured resources and longer-horizon scenario

| Quantity | Result | Interpretation |
| --- | ---: | --- |
| Successful selected attempts | 1,618 | Integrity checks passed; posterior not qualified |
| Unresolved selected attempts | 2 | Memory needs and sampler behavior unresolved |
| Successful worker wall duration | 1,612.11 hours | Sum of observed elapsed attempt durations, not CPU time |
| Successful saved-alignment size | 42.40 GiB | Actual file sizes; contents not freshly rehashed by this census |
| Proposed fresh chains | 1,620 | Complete original grid; different deterministic seeds |
| Proposed run length | 10,000 iterations | Unlaunched scenario, not a stopping rule |
| Linear successful duration scenario | 16,121.14 worker hours | Assumes unchanged per-iteration cost; excludes unknown failed outcomes |
| Linear successful alignment size scenario | 420.26 GiB | Uses 1,001 versus 101 saved alignments; excludes failed outcomes and other artifacts |

The proposed records preserve every original model identity, prior label,
configuration alias and chain role. Seeds derive from a versioned namespace,
horizon, chain identity and collision salt; all 1,620 are distinct and disjoint
from the original seed set. Each fresh attempt would retain 1,001 alignments
at ten-iteration intervals. Prespecified cutoffs of 2,500 and 5,000 would leave
750 and 500 retained alignment samples per chain, respectively. No original
sample concatenation, shortened failure horizon or tolerance relaxation is
proposed.

Neither linear scenario forecasts completion or convergence. The two failed
chains have no successful runtime baseline or qualified memory requirement.
Worker elapsed time includes the original scheduling and contention conditions.
Alignment complexity and posterior exploration can change over a longer run;
new native and downstream outputs would require additional storage. No native
memory peak is inferred from a receipt or an address-space cap.

The largest successful native alignment file is 1,392,681,456 bytes. Under the
same sample-size-only scenario it would reach about 13.8 GB, exceeding the old
two-GiB per-file cap. A future execution plan must revise output limits and
budgets explicitly. All proposed records currently have
`production_launch_allowed: false`.

BAli-Phy's [version 4.3 guide](https://www.bali-phy.org/README.html)
describes data-dependent convergence and recommends checking independent runs.
A longer iteration count cannot establish an adequate ancestral ensemble.
The closed scalar screens currently accept none of the complete quartets on
every scalar or candidate-length requirement.

## Numerical behavior across the full grid

The survey independently reads every logged iteration for five quantities:
RS07 rate, mean indel length, gamma shape, alignment width and total indel
length. Finite minima/maxima, original first/last tokens and every nonfinite
or out-of-support observation remain serialized. Zero total indel length is
allowed; other surveyed quantities must be positive. This is a numerical
survey, not a full prior-support or model-adequacy test.

Of the 1,618 selected intact chains, 134 contain a logged infinite gamma-shape
value during initialization: 54 broad, 27 centered and 53 package-prior chains.
Across all selected attempts there are 336 such observations, occurring only
at iterations 1–14. None occurs after either existing burn-in cutoff. These
events remain review evidence; they do not alone prove that the retained
posterior is corrupted, stationary or adequate.

Four intact selected chains record allocation warnings: one broad, two centered
and one package-prior chain. Across the 1,623 initial/selected-recovery attempts,
seven attempts have allocation warnings and five end with `std::bad_alloc`.
Those five failed attempts correspond to the three original failures and two
failed selected recoveries. Earlier additional historical attempts outside
this initial-plus-selected scope remain preserved.

Both unresolved selected chains belong to OG0000972. The domain/package chain
logs alignment widths 1,429, 46,264 and 376,981 at iterations 0, 1 and 2; its
last gamma shape is infinite. The whole/broad chain logs widths 10,873 and
491,177 at iterations 0 and 1, alongside rate 3.94e17 and mean indel length
about 55,779. Both fail even with a 48-GiB address-space cap. Their large
alignment expansions co-occur with extreme parameters; this does not isolate
causality or show that another memory cap would suffice.

The [four-panel figure](figures/baliphy_horizon_resources_20261003.pdf) shows
measured worker duration for all 13 families, distinct completed-chain reviews
by prior, every selected nonfinite event and both failed partial trajectories.
The [complete table](tables/baliphy_horizon_resources_20261003.tsv) retains all
405 quartets and their original chain IDs. Settings and chains share inputs;
these descriptive counts carry no independent biological replication claim.
PNG, SVG and PDF exports were inspected visually; source and artifact hashes
are in the [figure receipt](../metadata/baliphy_horizon_resource_figure_20261003.json).

## Reproducibility and execution limits

The [source survey](../scripts/baliphy_horizon_resource_inventory.py) freshly
hashes all used receipts, configurations, scalar logs, stderr and model/input
files. Closed recovery/scalar archive hashes bind the source lineage. Large
native alignment files are sized without newly hashing their contents; full
native residue/count verification is a separate workflow. Every original
attempt and all failed partial rows remain intact.

The full synthetic software workflow tested all 1,620 entries, 1,623 attempts,
405 quartets and fresh seeds. It exercised actual producer/readback commands
and refused completed restarts. Ten false serialized exports and two malformed
iteration traces were rejected. Synthetic journal fields are not evidence
of real process completion or a biological pilot. The final stable validation
is [version 2](../metadata/baliphy_horizon_resource_grid_validation_20261003_v2.json).

The prelaunch inventory covered 3,246 scalar/stderr files totaling 1.16 GB;
the largest is 740,649 bytes. Installed limits were two CPU equivalents,
eight GiB memory, no swap and one BLAS thread. Workspace estimate four GiB,
output allowance half a GiB and minimum free disk 64 GiB are planning values.
The 0.02–2-hour range per stage was uncalibrated planning, not an ETA.
Actual producer and reader used about 44 CPU-seconds each; controller journal
memory peaks do not establish native sampler memory requirements.

```bash
python scripts/check_baliphy_horizon_resources.py --output NEW_VALIDATION.json
python scripts/plot_baliphy_horizon_resources.py \
  --prefix NEW_FIGURE_PREFIX --table NEW_TABLE.tsv --receipt NEW_RECEIPT.json
```

Do not relaunch the completed census invocations. The next sampling design
requires sampler/initialization review, defensible memory/output budgets,
complete downstream diagnostics and the existing biological root/model/predictor
controls. The 10,000-iteration design remains unlaunched. Original full
native-residue and categorical readbacks continue unchanged; GPU inference
stays paused. All eight biological aims remain incomplete.


## Initialization follow-up and supported gamma limit

The pinned native source explicitly handles positive infinite gamma shape as
the equal-rate limit. The 336 initialization observations are retained as review
evidence; their values alone do not show invalid retained samples. The frozen
census and launch-disabled longer-horizon proposal are unchanged. A separate
[full reference-startup preflight](baliphy-reference-initialization-20261003.md)
now tests a supplied-alignment initial state using the original density and
alignment moves. This is startup verification, not a new posterior horizon or
a demonstrated repair of either earlier allocation failure.
