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
| 5. Ecological/morphological transitions | Evidence-curated traits and coverage diagnostics | Establish independently replicated usable transitions, preserve ambiguous assignments, and test controlled associations with structural change. A trait's number of labeled tips is not its number of independent origins. See [ecology](ecology-evidence-workflow.md). |
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
  completed: 1,000 fits total. Optimization/model review, exposure integration
  and conditional coupling remain downstream; see [coupling](conditional-site-coupling.md).
- Whole-protein duplicate/reference alignments are running with numerical
  auditors and order summaries queued behind them.
- Domain comparisons, common-core numerical verification, alternative-setting
  robustness, cross-guide comparisons and candidate sampling summaries completed.
  See [candidate sampling](duplication-domain-candidate-sampling-20260927.md).
  Matched backgrounds and phylogenetic association tests remain outstanding;
  none of these computational stages establishes a duplication effect.
- Species-tree sensitivity and database retrieval continue independently.
  GPU prediction remains paused under the user's current authorization.

Completion of a service is insufficient by itself: its source-bound receipt,
full intended output scope and relevant artifact/numeric checks must pass.
Retain failed/unestimable comparisons and denominator changes in downstream
analyses; do not silently shrink the scientific scope to successful cases.

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
