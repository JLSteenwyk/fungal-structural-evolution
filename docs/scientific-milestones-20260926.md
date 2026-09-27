# Open scientific milestones — September 26, 2026

This tracker preserves the eight aims in the [original objective](objective.txt).
**None of the eight aims is complete.** Completed computational stages below
support the aims but do not replace their statistical or biological requirements.
Chronological receipts and process records remain in [progress](progress.md).

| Aim | Current evidence | Next required result and completion evidence |
| --- | --- | --- |
| 1. Branch/clade structural change | Audited paired branch point estimates for 122 ESMFold and 125 expanded AlphaFold markers; direct-geometry/tree-path benchmark | Finish cohort-specific resampling and audits; integrate direct-coordinate and structural-alphabet evidence, model/coverage controls and shared ancestry; report calibrated branch/clade tests with multiple-testing treatment and uncertainty. See [recovery](recovery-20260926.md) and [geometry](tree-path-geometry.md). |
| 2. Sequence–structure coupling | Full ESMFold conditional site-rate models, whole-marker resampling and copy-omission sensitivity completed | Extend the qualified AlphaFold analysis; assess phylogenetic, model, alignment and prediction dependence; identify excess/reduced structural change conditional on sequence divergence. Fixed-rate conditional associations alone do not finish this aim. See [coupling](conditional-site-coupling.md). |
| 3. Domain evolution | Whole-proteome annotation/structural registry, clusters and family joins; duplicate/reference architecture controls | Complete within-domain comparisons; reconcile supported gains, losses, fusions, duplications and rearrangements separately from family turnover; propagate annotation and guide uncertainty. Annotation differences alone are not events. See [registry](whole-proteome-structure-domain-registry-20260923.md) and [controls](duplication-domain-controls-20260926.md). |
| 4. Duplication-associated change | Full terminal candidate/reference ledgers, guide sensitivity, sequence covariates and domain controls | Complete coordinate checks and all planned comparisons; join successful comparisons with explicit exclusions, examine reference/domain/confidence sensitivity, construct nonduplication backgrounds and account for family/phylogenetic dependence. Test divergence/asymmetry and functional-site changes. See [duplication](duplication-sister-references-20260926.md). |
| 5. Ecological/morphological transitions | Evidence-curated traits and independently verified ESMFold/expanded AlphaFold shared-site coverage | Establish independently replicated usable transitions, preserve ambiguous assignments, and test controlled associations with structural change. A trait's number of labeled tips is not its number of independent origins. See [ecology](ecology-evidence-workflow.md). |
| 6. Functional locations | Audited residue accessibility and functional correspondences | Join changes to supported core/surface and catalytic/binding annotations; evaluate matched backgrounds and supported interfaces/pockets where evidence permits. Annotation correspondence does not establish activity. See [functional workflow](functional-site-workflow.md). |
| 7. Selection | Codon fits and optimization/eligibility diagnostics | Resolve copy, alignment, saturation and optimization concerns; define justified test sets and multiple-testing scope; map supported residues. Structural acceleration is not evidence of positive selection. See [codon workflow](codon-model-environment.md). |
| 8. Ancestral/mechanistic cases | Exploratory families and functional candidates identified | Select cases supported by the preceding analyses, reconstruct ancestral amino-acid distributions and topology sensitivity, then predict authorized ancestral alternatives and formulate experimentally testable hypotheses. No completed ancestral reconstruction is claimed. |

## Current execution dependencies (updated September 27)

- Expanded AlphaFold resampling remains running across all 125 markers and
  75,000 paired draws, followed by automatic native-output/interval audit.
- ESMFold resampling and its full audit completed: 122 markers, 73,200 draw
  attempts, 73,167 estimable draws and 33 explicitly unestimable draws. These
  fixed-topology conditional intervals cannot substitute for AlphaFold intervals.
- Recovered AlphaFold Gamma4 and FreeRate4 site-rate exports and output audits
  completed: 1,000 fits total. Rate optimization/model review remains running. Accessibility, site summaries
  and functional annotation integration are complete; rate/exposure integration
  and conditional coupling remain downstream; see [coupling](conditional-site-coupling.md).
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
  complete. The full 144,040-fit phylogenetic working-model grid is running;
  a separate 432,120-fit ordinary-ML linear/quadratic/cubic grid is also running
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
  followed by full audit. Optimizer behavior is not a test of selection.
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
