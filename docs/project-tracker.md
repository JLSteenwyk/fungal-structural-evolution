# Project tracker

This tracker maps the full project objective to its current authoritative
evidence. “Prepared” and “running” are not completed scientific analyses.
None of the eight evolutionary aims is complete.

## Foundation and data custody

| Requirement | Current evidence | Status and next gate |
| --- | --- | --- |
| 501 fungal entries plus 25 outgroups | [Sampling manifest](../metadata/sampling_manifest.tsv), [analysis manifest](../metadata/analysis_manifest.tsv), [composition audit](../metadata/analysis_manifest_audit_20261006_v2.json) | The frozen analysis manifest has exactly 501 fungal entries and 25 outgroups across 19 explicit fungal phyla. One additional fungal candidate is retained only in the candidate manifest with its documented no-proteome exclusion; taxon quality, nomenclature and ecology remain explicit review variables. |
| Assemblies, annotations, proteins and coding sequence | [Original CDS audit](full-unmodified-cds-translation-20261006.md), [coding–structure source join](full-coding-structure-source-coupling-20261005.md) | Complete source custody and readback; annotation/code/compartment review still gates codon analyses. |
| Existing predicted structures | [Prediction-atlas union](full-prediction-atlas-union-20261005.md) | 2,961,055 source models inventoried; models are not automatically accepted for every downstream analysis. |
| Protein/domain structural-search infrastructure | [Full database workflow](full-atlas-foldseek-database-20261006.md) | Running original-CIF/native-coordinate verification; independent reader is queued behind it. |
| PAE uncertainty data | [Full PAE workflow](full-atlas-pae-20261005.md) | Full missing-AFDB retrieval is running; independent whole-queue reader is queued behind it. |
| Domain intervals and atoms | [Completed domain extraction](completed-domain-extraction-20261005.md), [atom readback](completed-domain-atom-readback-20261006.md) | Full export and retained-output integrity review are complete; biological domain boundary, confidence and homology checks remain. |

## Comparative framework

| Requirement | Current evidence | Status and next gate |
| --- | --- | --- |
| Species phylogeny, support and sensitivity | [Phylogenetic workflow](phylogenetic-workflow.md), [native taxon PMSF plan](../metadata/species_taxon_pmsf_plan_20261001.json) | Sixteen crossed taxon/alignment/guide PMSF refits run serially; do not select one guide before complete collection, comparison and support review. |
| Gene families, domains and reconciliation | [Orthology workflow](orthology-workflow.md), [domain workflow](domain-annotation-workflow.md) | Full family-discovery/reconciliation stages and independent readbacks remain required; structural clusters are not orthogroups. |
| Hydrophobin candidate inventory | [Provisional profile inventory](hydrophobin-profile-inventory-20261006.md), [immutable receipt](../metadata/hydrophobin_profile_inventory_20261006_v1.json) | All three named Pfam profiles were retained as candidates with taxon/protein links and source-model availability. Homology, reconciliation, coordinate qualification and structural/evolutionary inference remain required. |
| Dating | [Dating workflow](dating-workflow.md) | Relative divergence can be reported where warranted; calibrated branch-time analyses require defensible calibration and uncertainty review. |
| Ecological transitions | [Ecology evidence workflow](ecology-evidence-workflow.md) | Existing species-level statements and trait candidates are incomplete; tests require expanded curation, supported independent transitions and trait-uncertainty sensitivity. |
| Predictor and experimental controls | [Predictor-coordinate controls](predictor-coordinate-and-phylogeny-controls-20261004.md), [experimental-structure workflow](experimental-structure-workflow.md) | Existing controls inform later robustness checks; they do not establish prediction independence or biological effects. |

## Eight evolutionary aims

| Aim | Completion evidence required | Current gate |
| --- | --- | --- |
| 1. Branch/clade structural change | Reconciled homologous families, supported genealogy alternatives, direct-geometry benchmark, uncertainty and multiple-testing results | Structural atlas/database completion; reconciled trees; calibrated branch model. |
| 2. Sequence–structure coupling | Joint branch model with family/shared-ancestry effects, prediction/topology/alignment sensitivity and calibrated inference | Matched preliminary estimates exist; broad phylogenetic and atlas integration remain incomplete. |
| 3. Domain architecture and turnover | Independently supported gains/losses/fusions/rearrangements on reconciled families, distinct from within-domain change | Full structural clustering/search and family reconciliation remain open. |
| 4. Duplication and structural divergence | Reconciled duplication events, duplicate-pair controls, asymmetry model and mapped functional-site tests | Discovery and measurement controls are preparatory; no calibrated duplication effect is accepted. |
| 5. Ecology/morphology associations | Curated species traits, replicated supported transitions, phylogenetic model, trait/missingness sensitivity and multiplicity control | Trait evidence and same-method structural coverage are insufficient for inference. |
| 6. Functional localization | Core/surface/site/pocket/interface annotations, mapped structural changes and localization uncertainty | Requires eligible structural comparisons and independently supported feature annotations. |
| 7. Selection hypotheses | Eligible codon alignments/divergence, genetic-code and copy controls, selection model diagnostics and structural mapping | Source and translation evidence is ready; family-level eligibility and inference remain open. |
| 8. Ancestral reconstruction/case studies | Supported focal families, ancestral sequence ensembles, structure predictions, propagated uncertainty and proposed experiments | Computational diagnostics do not establish adequate ancestral posterior inference; focal cases remain to be selected after discovery safeguards. |

## Decision and completion rules

- The [research plan](research-plan.md) defines hypotheses, analysis order and robustness requirements.
- The [decision log and unresolved questions](decisions.md) records choices that must not be silently converted into biological conclusions.
- The [annotated bibliography](bibliography.md) records the literature basis for methods and focal biological controls.
- Completion requires evidence for every applicable row above, finalized reproducible outputs, figures/tables/methods, and mechanistic case-study proposals. An unavailable analysis needs an explicit, data-supported applicability determination; it is not silently treated as complete.
