# Open scientific milestones — updated September 28, 2026

This tracker preserves the eight aims in the [original objective](objective.txt).
**None of the eight aims is complete.** Completed computational stages below
support the aims but do not replace their statistical or biological requirements.
Chronological receipts and process records remain in [progress](progress.md).

| Aim | Current evidence | Next required result and completion evidence |
| --- | --- | --- |
| 1. Branch/clade structural change | Audited paired branch point estimates for 122 ESMFold and 125 expanded AlphaFold markers; direct-geometry/tree-path benchmark | Finish cohort-specific resampling and audits; integrate direct-coordinate and structural-alphabet evidence, model/coverage controls and shared ancestry; report calibrated branch/clade tests with multiple-testing treatment and uncertainty. See [recovery](recovery-20260926.md) and [geometry](tree-path-geometry.md). |
| 2. Sequence–structure coupling | ESMFold and recovered AlphaFold conditional site-rate models, whole-marker resampling and copy-omission sensitivity completed | Integrate source-specific results; assess phylogenetic, model, alignment and prediction dependence; identify excess/reduced structural change conditional on sequence divergence. Fixed-rate conditional associations alone do not finish this aim. See [ESMFold](conditional-site-coupling.md) and [AlphaFold coupling](recovered-afdb-coupling-20260928.md). |
| 3. Domain evolution | Whole-proteome annotation/structural registry, clusters and family joins; duplicate/reference architecture controls | Complete within-domain comparisons; reconcile supported gains, losses, fusions, duplications and rearrangements separately from family turnover; propagate annotation and guide uncertainty. Annotation differences alone are not events. See [registry](whole-proteome-structure-domain-registry-20260923.md) and [controls](duplication-domain-controls-20260926.md). |
| 4. Duplication-associated change | Full terminal candidate/reference ledgers, guide sensitivity, sequence covariates and domain controls | Complete coordinate checks and all planned comparisons; join successful comparisons with explicit exclusions, examine reference/domain/confidence sensitivity, construct nonduplication backgrounds and account for family/phylogenetic dependence. Test divergence/asymmetry and functional-site changes. See [duplication](duplication-sister-references-20260926.md). |
| 5. Ecological/morphological transitions | Evidence-curated traits and independently verified ESMFold/expanded AlphaFold shared-site coverage | Establish independently replicated usable transitions, preserve ambiguous assignments, and test controlled associations with structural change. A trait's number of labeled tips is not its number of independent origins. See [ecology](ecology-evidence-workflow.md). |
| 6. Functional locations | Audited residue accessibility and functional correspondences | Join changes to supported core/surface and catalytic/binding annotations; evaluate matched backgrounds and supported interfaces/pockets where evidence permits. Annotation correspondence does not establish activity. See [functional workflow](functional-site-workflow.md). |
| 7. Selection | Codon fits and optimization/eligibility diagnostics | Resolve copy, alignment, saturation and optimization concerns; define justified test sets and multiple-testing scope; map supported residues. Structural acceleration is not evidence of positive selection. See [codon workflow](codon-model-environment.md). |
| 8. Ancestral/mechanistic cases | Thirteen exploratory families; all 8,708,760 refined whole-protein and 3,365,640 refined domain amino-acid probabilities numerically checked. Refined/baseline, bound and full alternate-start comparisons completed; 26,126,280 alternate probabilities independently checked. | Resolve joint insertion/deletion and alignment uncertainty, posterior convergence, model adequacy and root/topology sensitivity; qualify cases using preceding biological analyses, predict authorized ancestral alternatives and formulate testable hypotheses. Conditional amino-acid marginals and 20-step diagnostic chains are not final ancestral ensembles. See [ancestral evidence](ancestral-case-inputs.md). |

## Current execution dependencies (updated September 28)

- Expanded AlphaFold resampling and native-output audit completed across all
  125 markers and 75,000 paired draws (23 unestimable). A separate
  [warning census and interval-table review](afdb-paired-resampling-review-20260928.md)
  completed; conditional intervals do not establish evolutionary significance.
- ESMFold resampling and its full audit completed: 122 markers, 73,200 draw
  attempts, 73,167 estimable draws and 33 explicitly unestimable draws. These
  fixed-topology conditional intervals cannot substitute for AlphaFold intervals.
- Recovered AlphaFold Gamma4 and FreeRate4 site-rate exports and output audits
  completed: 1,000 fits total. All 2,000 optimization diagnostics, selected
  rate-model comparison and full readback completed. The rate/exposure frame
  covers 125 markers and 47,529 sites. Full and copy-omission coupling fits
  and resampling completed, with exact subset and all 72 focal omission
  coefficients checked against leave-one-marker-out estimates. All 24
  specifications and both cohorts appear in the verified
  [coupling figure and interpretation](recovered-afdb-coupling-20260928.md).
- FastML independent replay covers all 153 nonempty inputs and 8,058,340
  probabilities. The targeted branch-scaling cache correction passed six
  paired regression runs; the full two-variant grid and separate numerical
  boundary-aware replay completed for 306 nonempty fits. Three empty inputs per variant remain
  explicit. Optimization and model qualification remain unresolved; see the
  [discrepancy record](fastml-indel-likelihood-discrepancy-20260927.md).
- Extant residue-anchor projection passed all 135 effective-input fixtures,
  covering 405 saved alignments and 1,620 reconstructed candidate sequences.
  The full BAli-Phy chain grid remains running. Scalar and length screens are
  queued; cross-chain ancestral-state and alignment mixing still require
  assessment. [Coordinate extraction](ancestral-residue-anchors-20260928.md)
  does not qualify posterior samples. The first seven production chains reached
  1,000 iterations; their 707 saved alignments now have verified geometry and
  [ancestral-state traces](ancestral-state-traces-20260928.md). No four-chain
  quartet is complete in that frozen subset, so convergence remains unproven.
  Automatic state extraction now covers the full 1,620-chain design; the first
  seven corrected handoffs exactly match the independently verified arrays.
  The [full categorical diagnostic queue](ancestral-categorical-diagnostics-20260928.md)
  is now active across all 405 four-chain models. At launch, 163 extracted
  chains were available but no quartet was complete; no convergence claim.
- Primary whole-protein production finished all 412,800 dispositions. The strict
  audit failed on short alignments; a full diagnostic isolated 27 RMSD discrepancies
  to two-residue mappings. All 327 short mappings passed separate analytic checks
  and remain excluded from unique-rotation interpretation. Geometry production
  and full independent readback passed all 387,646 successful alignments. The
  327 degenerate mappings exactly match the short-alignment census. Usable-order
  production and full readback passed all 206,400 pair/mask rows. Full residue-
  correspondence sensitivity and a [descriptive figure](whole-protein-order-sensitivity-20260927.md)
  are independently checked across 193,642 primary pair/mask comparisons. The
  failed strict audit remains preserved. Reference geometry/order sensitivity is
  verified. Whole-protein common-core fits passed independent checks for all
  563,808 records; all 36,944 event/reference links now retain coverage eligibility
  across masks, orders and mapping definitions. See [triplet coverage](whole-protein-common-core-comparisons-20260927.md).
  Background alignments and biological duplication tests remain open.
- Domain comparisons, common-core numerical verification, alternative-setting
  robustness, cross-guide comparisons and candidate sampling summaries completed.
  See [candidate sampling](duplication-domain-candidate-sampling-20260927.md).
  Metadata-matched backgrounds, balance, selection shifts and taxon/family
  concentration are now independently verified across all 54 sensitivity
  scenarios; see [matching diagnostics](duplication-control-balance-20260927.md).
  Background domain measurements and all 82,944 descriptive summary settings are
  complete. The full 144,040-fit phylogenetic working-model grid, analytic-gradient
  checks and export completed; 658 unique fits remain for numerical review.
  A separate 432,120-fit ordinary-ML linear/quadratic/cubic grid is running
  on verified identical observations across five trees. The expanded designs and
  input inventory passed full readback; nonlinear joint support is verified for
  all 57,616 inputs, including separate nonnegative certificates for 112 edge cases. Full likelihood comparison exports retain all settings and fit flags.
  Numerical refinement, uncertainty and biological association tests remain outstanding;
  none of these computational stages establishes a duplication effect.
- Recovered AlphaFold accessibility and full site-summary readbacks are complete:
  30,618 models, 9,453,757 paired observations and 47,529 retained sites. Functional
  annotation integration preserves every site; 6,444 annotation rows map to 50
  paired sites, with conserved candidates at 18. The predictor comparison covers
  all 150 jointly observed exact protein coordinates and all 600 specified local
  coordinate comparisons. These inputs do not complete the functional aim.
- Local codon branch-profile checks cover all 1,632 cases. All eight saved starts
  per case (13,056 starts) are being reoptimized without the profile constraint,
  followed by full audit. A frozen partial review checks 777 completed cases:
  14 have between-start log-likelihood spreads greater than 10. The full source-
  model/parameter audit and 855 cases remain pending in that snapshot. Optimizer
  behavior is not a test of selection; see the
  [multiple-start review](local-mg94-unconstrained-multistarts-20260927.md).
- The September 27 cache refresh leaves 3,626 unique marker proteins without
  cached models. GPU prediction remains paused; this is a marker-gap count,
  not whole-proteome coverage or proof of database absence.
- Species-tree sensitivity and database retrieval continue independently.
  GPU prediction remains paused under the user's current authorization.

Completion of a service is insufficient by itself: its source-bound receipt,
full intended output scope and relevant artifact/numeric checks must pass.
Retain failed/unestimable comparisons and denominator changes in downstream
analyses; do not silently shrink the scientific scope to successful cases.

## Additional September 27 evidence

- Full joint-covariate support is verified for all 28,808 unique comparative
  inputs and 82,944 settings, including explicit nonnegative certificates for
  24 numerical edge cases. Conditional model fitting and inferential validation
  remain pending. See [joint support](joint-covariate-support-20260927.md).
- Exact marker/column checks found no cross-predictor overlap at the three
  currently mapped focal ecological taxa. Taxon-level coverage in both sources
  did not imply comparable marker membership. See [ecological edge coverage](ecology-optimal-edge-states-20260927.md).
- Eight taxa now have separately reviewed decay classifications and fully checked
  structural coverage. The primary two-tree mapping requires two undirected
  changes but no particular edge; uncertain coding changes localization. The
  full bootstrap mapping and independent checks completed across 2,000 trees and
  five codings. Primary unknown coding has no required-change edge on any tree;
  apparent localization depends on assigning Jaapia brown. See
  [decay uncertainty](wood-decay-phylogenetic-diagnostic-20260927.md).
- A [predictor geometry figure](figures/functional_predictor_geometry_20260927.pdf)
  retains all 600 local comparisons at 150 shared functional positions. It is
  a prediction-source control, not a map of evolutionary changes.

## Ancestral completion requirements (updated September 27)

- Whole-protein refinements and independent posterior readback are complete.
  The completed baseline/refined comparison changes the most probable amino
  acid in 14 of 435,438 dependent node/site comparisons, with no opposing
  calls supported at least 0.9 in both fits. Bound comparisons have no state
  changes. See the [figure](figures/whole_refinement_sensitivity_20260927.pdf).
- All 468 alternate-start whole-protein fits passed likelihood and probability
  checks. Full comparison is complete: 128 most-probable state changes across
  1,306,314 dependent node/site comparisons, none opposing at 0.9 support in
  both fits. The largest support shift occurs in a slightly poorer fit and
  remains documented; best-fit agreement alone does not establish robustness.
- Independent gap-character estimates violated mutual-exclusion constraints;
  conditioning them on compatibility is not a fitted joint evolutionary model.
  Historian's completed largest-family MAFFT retries also show minimum-edge
  sensitivity (10–15 ungapped edits at non-root candidates). Retain these
  limitations when evaluating candidate sequences.
- BAli-Phy short-chain output checks are complete for all 324 configurations
  and 3,888 candidate-node samples. All 405 explicit-prior initializations also
  passed. The full 1,620-chain grid (135 effective inputs × three priors × four
  seeds) is running with 16 CPU workers and a 1,000-iteration first horizon.
  All 405 quartet scalar checks are queued at both burn-in cutoffs. The complete
  short-run resource audit and isolated crash-recovery tests passed; initial
  long-chain workers and enforced limits were verified. Posterior mixing,
  alignment/ancestral-state diagnostics, model/root sensitivity and qualified
  ancestral samples are still required. Exact-input aliases retain their
  labels but never count as independent chains. See the
  [BAli-Phy assessment](baliphy-method-assessment-20260927.md).
- GPU prediction remains paused. No qualified ancestral structural ensemble
  has been produced, and aim 8 remains open alongside the other seven aims.

## Cross-cutting requirements

Uneven model coverage is substantial: only about 23.4% of the examined terminal
singleton-side duplication events have both models in the frozen AlphaFold
inventory, from 153 of 526 reconciled taxa. The 527-entry candidate manifest
includes excluded *S. jurei*. Use analyzed membership, source-specific model
availability and uncertainty explicitly; do not interpret missing models as
biological absence. See [coverage](duplication-sampling-coverage-20260926.md).

Branch lengths and structural distances have different meanings. Without
supported absolute calibration, report relative divergence rather than change
per unit time. [Dating evidence](dating-workflow.md), prediction-source
benchmarking, domain-orientation uncertainty, gene-copy review and multiple
testing remain requirements across the relevant aims.

The final dataset/atlas, phylogenies, uncertainty estimates, evolutionary tests,
figures, methods and mechanistic case studies must all retain their actual
sampling and validation scope. No completion date or overall percentage follows
from the fraction of computational jobs finished.

## September 27 experimental case-control update

All 706 nominated experimental coordinate entries passed full CA mapping/readback;
all 65,556 observed-coverage interval rows and 78 case summaries passed independent
checks. Candidate count reductions are substantial for some cases. Direct
four-protein fits, contrast readback and experimental-reference sensitivity
are now complete: 115,200 reference-unit records and 936 zero-inclusive case
summaries were checked. All 13 cases, four metrics, six coverage screens and
three margins are retained. The combined experimental/ancestral case table
preserves all 936 settings, including 648 with no qualifying reference units.
These are dependent descriptive comparisons; mechanistic interpretation and
the broader evolutionary aims remain pending. See
[experimental case controls](case-independent-control-coverage-20260927.md).

A literature-supported potential Neocallimastix species complex also motivates
six new representative-sensitivity guide searches, now running. This does not
establish a revised unique-species count or resolve the other uncertain labels.
See [taxon identity](taxon-identity-sensitivities.md).

See the [restored-access checkpoint](restored-access-checkpoint-20260928.md)
for terminal-state and artifact verification of the newly completed stages.

September 28: [analytic refinement](matched-reml-analytic-refinement-20260928.md)
is active for all 658 remaining flagged working-model fits. Dense-likelihood
and derivative-free synthetic checks passed; production qualification remains
pending output readback and review.

September 28: the [full paired FastML comparison](fastml-paired-comparison-20260928.md)
checked all 8,058,340 paired marginal entries. Cache refresh changes eight
OG0000294 fits and can worsen the attained likelihood. Native gamma rates
differ from SciPy rates at all 306 fitted shapes. The full native-rate replay
and source audit completed: maximum likelihood error in refreshed fits is
1.82e-8; the three uncorrected stale-cache errors remain. The two worse optima
still require optimization review.

September 28 optimizer audit: all 306 nonempty FastML variant fits used an
effective one-model-iteration limit and emitted the corresponding limit
message. Full five-start refinement inputs are prepared for all156 designs
(765 nonempty fits, 15 empty dispositions). Execution is active in version 3,
with explicit protection against native tree annotation and branch-floor
behavior. Earlier failed handoffs are preserved.
See the [optimizer evidence](fastml-paired-comparison-20260928.md).

September 28: full658 matched-model refinement readback is queued after
producer completion. The checker passed a frozen152 completed records plus
three tamper cases. Whole-grid qualification and setting-map integration
remain pending.
