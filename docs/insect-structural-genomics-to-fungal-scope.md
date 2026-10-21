# From comparative insect structural genomics to fungal structural evolution

Wu et al. (2026) provide a useful precedent for combining a broad structural
atlas with phylogenetic context, remote-homology discovery and functional
annotation. The study motivates the atlas and discovery components of this
project; it does not supply fungal taxa, fungal family histories, a species
tree, or evidence for fungal evolutionary effects. The annotated citation is in
[the bibliography](bibliography.md).

## What transfers directly

| Precedent component | Fungal implementation | Required evidence before interpretation |
| --- | --- | --- |
| Broad comparative structural catalogue | Trace all source AFDB/ESMFold models, including alternative products and missingness | [Prediction-atlas union](full-prediction-atlas-union-20261005.md), coordinate/PAE qualification and source-specific controls |
| Structural search and candidate remote homology | Full-atlas Foldseek database and later candidate partition, followed by database/function review | [Clustering handoff](full-atlas-structural-clustering-handoff-20261006.md), direct alignment checks and family/reconciliation review |
| Functional annotation of structurally related proteins | Sequence, domain and structural annotations remain separate evidence streams | Supported annotation evidence; no transfer of experimental function from structural similarity alone |
| Phylogenetic context | A 501-fungus, 25-outgroup sampling frame with a sequence-based species-tree sensitivity collection | [Sampling record](outgroup-selection-20261006.md), completed PMSF collection, discordance and rooting review |

## Fungal extensions

The project adds the following analyses because a structural atlas alone cannot
identify the evolutionary process behind an apparent fold or cluster pattern.

1. **Branch-specific sequence–structure coupling.** On supported homologous
   families and genealogy alternatives, estimate amino-acid and 3Di changes,
   then test structural excess conditional on sequence change, branch duration,
   alignment coverage, protein length, family effects and prediction uncertainty.
   Structural-alphabet branch lengths are not interpreted as physical
   displacement; direct local/contact geometry is the benchmark.

2. **Reconciliation-aware duplication and domain evolution.** Analyze gains,
   losses, fusions, rearrangements and duplicate asymmetry separately from
   within-domain structural change. Structural clusters remain candidate
   similarity groups, not orthogroups or duplication events.

3. **Genome architecture integration.** Combine qualified structural/family
   results with retained gene-order, repeat, collinearity and genome-context
   covariates only after RepeatModeler/RepeatMasker and reconciliation gates
   complete. This permits tests of whether architecture changes and sequence or
   structural changes co-occur on supported branches without treating genomic
   proximity as causal evidence.

4. **Ecological and morphological transitions.** Test only curated traits with
   repeated, independently supported changes on an accepted phylogeny, while
   propagating uncertain or mixed assignments and correcting for multiplicity.

5. **Functional localization and selected case studies.** Map accepted changes
   to cores, surfaces, pockets, catalytic residues and supported interfaces.
   Promote families to ancestral reconstruction only through the predeclared
   [case-study protocol](case-study-prioritization-20261006.md).

## Interpretation safeguards

Predicted structure derives from sequence and therefore is not an independent
measurement of evolution. Leading findings must be evaluated across prediction
sources/settings, confidence and PAE masks, topology and alignment alternatives,
and experimental structures where available. Low-confidence/disordered proteins
and uncertain domain orientations remain explicit dispositions rather than
being silently removed or called novel folds.

The operational completion gates and the status of all eight evolutionary aims
are maintained in the [project tracker](project-tracker.md). The project
[research plan](research-plan.md) specifies the analysis sequence, uncertainty
requirements and completion criteria. This note is a scope translation and not
a completed analysis or a claim about fungal structural evolution.
