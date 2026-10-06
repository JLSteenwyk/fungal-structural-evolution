# Mechanistic case-study prioritization — 2026-10-06

This protocol controls promotion of a protein family or domain from discovery
output to a mechanistic case study. It applies to the 13 existing exploratory
ancestral-case families and to any new candidate emerging from the full atlas.
It does not select a case now. No candidate has yet met these requirements.

## Eligibility gates

A candidate must pass every applicable gate before it can be ranked for a
case study.

| Gate | Required evidence | Disqualifying shortcut |
| --- | --- | --- |
| Family identity | Homology assessment, reconciled gene tree, duplicate/loss uncertainty, and retained taxon/protein identities | Calling a Foldseek cluster, Pfam hit, or structural group an orthogroup |
| Species-tree robustness | Result persists across supported taxon/marker/model sensitivity views, with affected branch/path mapping reported | Assigning a merged or root-dependent path to one branch |
| Structural change | Direct geometry and structural-alphabet estimates are concordant in direction or explicitly explain their disagreement; domain orientation is separated from within-domain change | Treating a structural-alphabet length as Å displacement |
| Prediction uncertainty | Coordinate confidence/PAE qualification, alternative-predictor sensitivity, and experimental-reference evidence where available | Treating high confidence or predictor agreement as experimental validation |
| Sequence context | Matched sequence divergence and alignment coverage are reported; a selection test is run only on an eligible coding alignment | Calling structural acceleration positive selection |
| Domain and copy context | Gains, losses, fusions, rearrangements, and duplication events are separately reconstructed from within-domain changes | Treating a changed domain architecture as a changed homologous fold |
| Quality and ecological context | Source/annotation/FCS/codon/identity flags are reviewed; ecological claims require curated, independently replicated transitions | Inferring an ecological mechanism from one species or an unreviewed label |

Missing evidence remains a visible disposition. A case may be retained as a
discovery candidate but cannot be ranked as a mechanistic result by substituting
an unavailable gate with a favorable proxy.

## Priority evidence package

After the gates are met, report an evidence table for every eligible candidate
instead of an opaque composite score. Each row must record: family/domain ID;
event and branch/path IDs; all sensitivity memberships; sample sizes; direct
geometry and structural-alphabet summaries; sequence divergence; predictor and
experimental-reference coverage; domain/copy state; quality flags; functional
annotation source; and explicit uncertainties. Rank candidates only within
their comparable analysis stratum, using these prespecified priorities:

1. Replicated structural signal across independent analytical sensitivities.
2. Clear separation of within-domain change from architecture or copy-number
   turnover.
3. A compact, well-supported clade that permits extant and ancestral sequence
   ensembles without discarding conflicting alternatives.
4. Supported functional annotation or experimentally tractable biochemical
   phenotype, stated at the level supported by its source.
5. A feasible validation experiment that tests a concrete prediction.

An apparent large displacement, unusual fold, extreme branch length, or sparse
lineage alone is not a priority criterion. The 13 existing candidate families
remain exploratory because reconciliation, accepted species framework,
full-domain/whole-protein comparison and ancestral posterior requirements are
still incomplete.

## Validation-experiment proposals

For each promoted case, formulate a falsifiable experiment only after the
evidence package identifies a supported function and affected feature. The
proposal must compare extant and uncertainty-sampled ancestral or reconstructed
variants, include expression/folding controls, and state the predicted result
that would falsify the model. Depending on the supported feature, appropriate
assays may include:

| Supported feature | Example testable proposal |
| --- | --- |
| Core packing or stability change | Compare purified variants with thermal/chemical unfolding and aggregation measurements. |
| Catalytic-site or substrate-pocket change | Measure activity and substrate specificity using a prespecified substrate panel and inactive-site controls. |
| Ligand/interface change | Measure binding or complex formation with reciprocal partner and concentration controls. |
| Domain fusion, loss, or rearrangement | Compare full-length and domain-resolved constructs for folding, localization, or regulated activity, with matched expression controls. |
| Surface or secretion-associated change | Test localization or interaction phenotype with signal-peptide/domain controls; do not infer phenotype from a surface map alone. |

These are future proposals, not completed experiments. Experimental structures,
when available, are a benchmark for coordinates and do not by themselves test
the proposed biological mechanism.

## Current evidence and next gate

The [ancestral-case inputs](ancestral-case-inputs.md) preserve 13 exploratory
families and their alignment/domain sensitivities. The [whole/domain inspection
cases](whole-domain-inspection-cases-20260927.md) retain descriptive contrasts
without function-based selection. The [hydrophobin profile inventory](hydrophobin-profile-inventory-20261006.md)
is a separate candidate profile census, not a case-study selection. Case
promotion awaits the supported species-tree collection, reconciled family
analyses, structural-atlas qualification, and ancestral ensemble inference.
