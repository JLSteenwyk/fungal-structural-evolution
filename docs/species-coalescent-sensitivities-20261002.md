# Full gene-tree species-tree sensitivities

Thirty ASTRAL-III 5.7.8 analyses and their complete independent local/global
quartet validation are complete. They use all 125 profile marker trees or all 125 MAFFT marker
trees, each under five actual taxon cohorts and three gene-support settings.
This estimates candidate species relationships from gene trees and complements
the original concatenated PMSF analyses. It does not accept a final species
tree, biological root or cause of discordance. October 2 artifact dates use UTC.

ASTRAL optimizes agreement with gene-tree quartets under a constrained search.
Its multispecies-coalescent and local support assumptions must be assessed
before biological interpretation. It accepts missing taxa and unresolved
gene-tree branches. References: [Zhang et al. 2018](https://doi.org/10.1186/s12859-018-2129-y),
[Sayyari and Mirarab 2016](https://doi.org/10.1093/molbev/msw079), and the
[official tutorial](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md).

## Complete data design

| Cohort | Taxon entries | Ingroup | Outgroups |
| --- | ---: | ---: | ---: |
| Full primary | 526 | 501 | 25 |
| Exclude sparse boundary taxon | 525 | 500 | 25 |
| Exclude taxa below 10% occupancy in both matrices | 522 | 499 | 23 |
| Exclude curated hybrids | 524 | 499 | 25 |
| Exclude hybrids and incomplete labels | 503 | 478 | 25 |

The four retained cohorts are the exact closed sets used for native PMSF
sensitivities. They are not redesigned after inspecting coalescent output.
The primary working set remains 526 entries/25 outgroups; 501 entries must
not be equated with 501 independently established species. Individual gene
trees need not include every outgroup. Full prelaunch source verification
passed 1,305 bindings, all 250 gene trees and every cohort union. Each retained
marker has 402–504 taxa depending on the cohort; no marker is discarded.

The three support policies are uncontracted point trees, contraction of
reported SH-aLRT below 10, and contraction of reported SH-aLRT below 80.
Values exactly at 10 or 80 are retained. Both filtered policies also contract
edges with unavailable support. All 565 profile and 366 MAFFT unreported
support values remain explicitly classified as missing, rather than zero.
The baseline preserves their original point-topology resolutions.

These are SH-aLRT sensitivity settings. The ASTRAL literature's low-bootstrap
contraction recommendation does not make SH-aLRT numerically interchangeable
with bootstrap support. No policy is preferred or treated as a posterior
gene-tree ensemble. Contraction occurs before cohort taxon pruning, so a
collapsed/pruned path does not acquire an invented inherited support.
Input branch lengths and labels are omitted because the native calculation
uses unrooted topologies. Original sequence-tree lengths remain in their
source archives.

## Input proof and native execution

There are 3,750 marker/cohort/support input states and 1,781,685 original
gene-split support decisions. The producer uses Bio.Phylo contraction and
pruning; a separate DendroPy reader checks every expected surviving split
by projecting the independently audited eligible source splits onto the
retained gene tips. It also checks every missing/low/exact-cutoff disposition,
marker order, input hash, tip union and ingroup/outgroup role. Full five-cohort/
three-policy synthetic tests passed and rejected 17 foreign-taxon/threshold
errors. These are software checks, not a biological pilot or production proof.

Input closure requires all 30 cases, full source/output hashes and both exact
original producer and reader completion journals. A new closure helper
handles legitimate null `MESSAGE` fields in structured journals; its regression
replayed all source hashes and the actual original completed gCF figure
journal. Existing launched helpers and their records remain unchanged.

The full input proof is complete: all 30 cases, 3,750 transformed gene trees,
1,781,685 source support decisions and 1,545,528 retained internal split states
passed the independent reader. Closure checked 1,408 bindings and both original
producer/reader completion journals. [Completed input locator](../metadata/species_coalescent_inputs_completed_20261002.json).

The installed ASTRAL jar and all three dependencies exactly match the official
5.7.8 archive. Native synthetic resolved, missing-taxon and polytomous inputs
passed with Java 21. Exact executable/library hashes, help/version output and
the archive hash are recorded outside Git and bound by the launch plans.
This checks software integrity and I/O; it does not validate the full native
species-tree estimates.

```bash
python scripts/prepare_species_coalescent_inputs.py \
  --plan metadata/species_coalescent_input_plan_20261002.json
python scripts/readback_species_coalescent_inputs.py \
  --plan metadata/species_coalescent_input_plan_20261002.json \
  --output results/phylogeny/full-species-coalescent-inputs-20261002-v1/readback.json
python scripts/close_species_coalescent_inputs_v2.py \
  --plan metadata/species_coalescent_input_completion_plan_20261002.json
python scripts/run_species_coalescent_sensitivities.py \
  --plan metadata/species_coalescent_native_plan_20261002.json
```

Commands identify the original frozen runs. Reproduction needs new plans and
new outputs; launched scripts/plans must not be changed. The completed small
input locator is `metadata/species_coalescent_inputs_completed_20261002.json`.
Native output is `results/phylogeny/full-species-coalescent-sensitivities-20261002-v1`.
All full source dictionaries, input trees, native trees and logs stay outside Git.

Every native case uses the pinned seed and `-t 2` to retain quartet fractions,
weighted quartet evidence, local posterior support, effective gene count and
internal branch lengths in coalescent units. Exact child PID/create/command,
return code and output hashes are recorded. Arbitrary native Newick rooting
does not assign biological polarity. Local posterior support is conditional
on the adjacent clades and gene/MSC model, not calibrated global topology
confidence. Effective gene counts can be fractional with unresolved inputs.
Coalescent lengths are not substitutions, elapsed time or structural displacement;
terminal lengths remain unavailable.

Resources were estimated before launch: each serial input stage uses two CPU
cores/16 GiB/no swap/one BLAS thread, with 2 GiB output allowance and 100 GiB
disk reserve. Each serial native case has two available processors, a 28 GiB
Java heap and 32 GiB service cap, no swap/GPU/new charges. ASTRAL-III is mostly
single threaded. The native 0.05–6 hours per case/1.5–180 hours per batch ranges
are uncalibrated planning estimates and exclude dependency waits; they are
not ETAs. Existing analyses continue and GPU prediction remains paused.

## Independent local and global quartet verification

All 30 native inferences and the corrected full readback are complete. It retains all 15,510
internal branches and all 1,938,750 branch–gene states across the 30 cases,
plus 3,750 per-gene global matching/resolved-quartet counts. It uses separate
Python/Numba algorithms and DendroPy parsing, without ASTRAL libraries.

For local evidence, a postorder dynamic program counts color-distinct quartets
by their leaf distribution among the children of their lowest common ancestor.
It distinguishes resolved 3+1, 2+2 and 2+1+1 arrangements from unresolved
1+1+1+1 arrangements. Resolved rooted triples retain their cherry pair for
the 3+1 case. Each integer count is divided by the number of available
four-clade combinations in that gene before pooling across genes. Missing
incident clades contribute unavailable evidence; unresolved quartets do not
become supporting observations. Effective gene counts can consequently be
fractional. The two NNI alternatives are compared as an unordered pair;
their original native labels and separate canonical group counts are retained.

Global matching quartets are calculated separately by summing compatible
node-component intersections over every species/gene internal-node pair.
Each matching resolved quartet has two endpoint centers in each tree and
is counted four times; division by four gives the exact numerator. An
independent denominator calculation subtracts unresolved quartets from all
four-tip combinations using the fourth elementary symmetric polynomial
of incident-component sizes at each multifurcation. Both methods preserve
missing taxa and polytomies. All arithmetic is integer before normalization.

Local posterior values are recomputed using beta-tail integration and log
normalization under the native default Yule prior, lambda = 0.5. Branch
lengths are recomputed as the corresponding MAP estimates in coalescent
units. Zero evidence retains unavailable quartet fractions and prior-only
posterior support. Native floating-point values are compared with absolute
and relative tolerance 2e-8. Numerical agreement does not establish model
adequacy, gene independence, biologically accepted rooting or search optimality.

Exhaustive four-tip enumeration passed 18,720 color/tree cases. Five actual
native contracts cover resolved, fractional, missing, zero-evidence and
discordant inputs; 120 altered native counts, posteriors or lengths were
rejected. Global scoring passed all 105 six-tip species topologies crossed
with 12 resolved/partial/multifurcating gene cases, plus four twelve-tip cases
(1,264 tree pairs). The initial software-only DendroPy API failure is preserved
in the v1 records; v2 specifies one input source and passed the full contracts.

The actual first full 526-taxon profile/uncontracted candidate passed all
523 branches and 65,375 local states. Independent global matching score
240,671,603,016 and resolved denominator 268,840,943,689 exactly match native
output (normalized score 0.8952193059343416). Its one original auditor completion
journal and 1,434 source/artifact bindings are closed in the
[named-case locator](../metadata/species_coalescent_quartet_preflight_completed_20261002_v2.json).
This is complete verification of one full candidate, not full-batch completion,
a biological pilot or a preference for this candidate.

The first full global/local readback took 33.2 seconds. Before queuing the
full serial readback, resources were set to two CPU/16 GiB/no swap/one BLAS
thread, a 4 GiB output allowance and 100 GiB disk reserve. Its 0.25–2 hour
planning range excludes native dependency waits. A separate closure follows
all 30 cases and requires both original native/auditor journals and full hashes.
No GPU or paid resources are used.

```bash
python scripts/readback_species_coalescent_quartets_v3.py \
  --plan metadata/species_coalescent_quartet_plan_20261002_v2.json
python scripts/close_species_coalescent_quartets.py \
  --plan metadata/species_coalescent_quartet_completion_plan_20261002_v2.json
```

These commands identify immutable completed runs, not instructions to start
duplicates. Full readback output is
`results/phylogeny/full-species-coalescent-quartet-readback-20261002-v2`;
the [completed small locator](../metadata/species_coalescent_quartets_completed_20261002_v2.json)
binds 1,757 source/artifact hashes and both original native/auditor completion
journals. All 15,510 branches, 1,938,750 local branch–gene states and 3,750
global per-gene states passed, including the previously failing case.
The [runtime checkpoint](../metadata/project_runtime_checkpoint_20261002_v6.json)
records 30/30 native cases and full numerical closure, 22 live pipeline handles,
six original scientific/retrieval jobs, 26 verified terminal successes and ten
preserved failures. It checked 1,874 distinct closed bindings and recorded
187,136/298,848 new structural comparison states at its observation time.
Later native progress remains outside Git.

## Remaining qualification

Input preparation, full independent input readback and closure are complete.
Native inference is complete with all source/output hashes and its original
completion journal verified; corrected full numerical/quartet readback and
two-journal closure are complete. Candidate comparisons against all concatenated
references are [complete across all 70 views/315 pairs](coalescent-reference-comparisons-20261002.md),
including independent raw-tree readback and inspected figures. Comparisons using
the 16 actual subset PMSF refits are [queued separately](native-subset-tree-comparisons-20261002.md).
Gene-tree
estimation uncertainty, locus dependence, low taxon occupancy, orthology,
gene-selection sensitivity, hybrid ancestry and MSC adequacy remain open.
These candidate topologies must be qualified before reconciliation and
structural-rate comparisons. All eight evolutionary aims remain incomplete.

## Native completion and effective-N correction, October 2 UTC

All 30 native cases completed with zero return codes. Every case/command,
per-case receipt and full source/output hash was checked, and the original
wrapper completion/resource journal was verified. The
[native output inventory closure](../metadata/species_coalescent_native_outputs_completed_20261002.json)
binds 1,571 hashes and one original journal. It certifies completed inference
outputs, not independently verified quartet numerics or biological eligibility.

The original full numerical reader passed the first case and stopped in the
second: native EN was 125.0, whereas independently pooled resolved evidence
was 124.99995265575164. Inspection of the
[author implementation](https://github.com/smirarab/ASTRAL/blob/master/main/phylonet/coalescent/WQInference.java)
shows that EN starts with the number of genes having all four incident clades
available. It is replaced by the fractional resolved total only when their
absolute difference exceeds 0.001. Quartet fractions, local posterior support
and MAP lengths use this native effective denominator. Raw fractional resolved
evidence and available-gene counts must therefore remain distinct.

The immutable v3 numeric module implements that rule and retains both counts
and the adjustment disposition for every branch. It uses the unchanged exact
quartet DP; comparison tolerances remain absolute/relative 2e-8. The installed
5.7.8 jar passed three actual fixed-reference contracts with exactly one
unresolved quartet among 999, 1000 and 1001 available combinations. These
cover both sides and the floating-point representation near the 0.001 rule.
All five prior native contracts (15 branches) passed again; all 120 changes
of 0.0001 to counts, effective N, posterior support or lengths were rejected.
Source/master inspection is supplemented by these actual installed-jar tests;
no package/tool update or inferred tree mutation was made.

The original failed full reader and its partial exports remain unchanged.
New v2 plans/output directories run the v3 reader over all 30 candidates.
The previously failing complete 526-taxon SH10 candidate has now passed all
523 branches and 65,375 local states, with its exact global matching score.
The complete production batch passed and its two-original-journal closure
binds all 1,757 hashes. All 315 corrected comparisons and their figures also
passed independent readback and closure; both actual figure pages were inspected.
The corrected comparisons and figure used new versioned dependencies. Their
handoff wrapper verifies original invocation-linked completion/resource
journals and rejects collected failed-unit default success fields. All eight
scientific aims, root/model/gene/taxon/dating/reconciliation qualification
remain open; GPU structure prediction stays paused.
