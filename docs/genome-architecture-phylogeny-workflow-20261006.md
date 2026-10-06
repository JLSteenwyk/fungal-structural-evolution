# Genome architecture and protein-structure evolution extension

This is a full-panel extension for testing whether changes in genome context
coincide with protein structural change. It does not report synteny, repeat,
rearrangement, or association results yet. Protein/domain architecture and
chromosomal genome architecture are distinct measurements and will remain so.

## Questions

1. Which branches show elevated disruption or conservation of orthology-anchored
   gene neighborhoods and syntenic blocks?
2. Do within-domain structural-change residuals covary with neighborhood
   disruption after sequence divergence, family, branch duration, prediction
   uncertainty, assembly contiguity and shared ancestry are accounted for?
3. Are apparent associations explained by repeat-rich context, gene duplication,
   annotation fragmentation, or family turnover rather than change inside a
   homologous protein domain?

This extends the project's domain-architecture aim; it does not redefine
structural change as a genome rearrangement.

## Inputs and preconditions

The analysis will use the complete 526-taxon genomic DNA and annotation
registries, assembly/FCS quality records, reconciled gene families, a supported
species-tree collection, and confidence-qualified structural measurements. The
current whole-genome and annotation source records preserve the necessary
coordinates, but source availability does not yet establish a synteny block or
repeat annotation. Fragmented assemblies, missing coordinates, unplaced
scaffolds, FCS-flagged loci, alternative products and uncertain taxon identity
remain explicit dispositions.

No family may receive a genome-context score before its gene/protein mapping
and reconciliation are qualified. No species-tree branch rate may be accepted
before the taxon/marker/model sensitivity collection is evaluated.

## Full-panel workflow

1. **Coordinate and assembly audit.** Build an immutable gene-order table from
   all retained annotation features, preserving scaffold, strand, ordinal
   position, intergenic distance, alternative-product and missing-coordinate
   status. Stratify every later result by assembly level/contiguity and FCS
   context rather than excluding difficult taxa by default.
2. **Repeat-context evidence.** Obtain repeat annotations with a versioned,
   taxon-consistent method and report library/model provenance, masked fraction,
   repeat class and distance-to-gene. A repeat annotation is context evidence,
   not proof that a rearrangement or gene birth was repeat-mediated.
3. **Orthology-anchored neighborhoods.** With reconciled families, enumerate
   local conserved-neighbor and collinear-block candidates. Retain one-to-many,
   unplaced and ambiguous-anchor dispositions; do not force a single anchor for
   duplicated genes or treat a shared Pfam domain as an orthology anchor.
4. **Phylogenetic mapping.** Model neighborhood/block disruption, local gene
   density, repeat proximity and architecture turnover on each supported
   species-tree alternative. Map uncertainty to paths when a gene tree cannot
   identify a unique species-tree branch.
5. **Joint structural model.** For homologous domains, model structural-change
   residuals given sequence change against genome-context variables using
   family and phylogenetic effects. Include branch duration when defensible,
   assembly contiguity, annotation density, protein/domain length, coverage,
   prediction source/confidence, duplication state and multiple-testing control.
6. **Sensitivity and interpretation.** Repeat supported analyses after
   excluding or down-weighting fragmented, FCS-flagged, repeat-ambiguous and
   taxonomically uncertain inputs. Test within-family and between-family
   effects separately. An association cannot establish repeat causality,
   adaptation, or a physical effect of a rearrangement on protein structure.

## Resource and launch gates

No full synteny or repeat discovery job is launched by this document. Before
any expensive stage, record the all-taxon input census, per-genome gene and
scaffold counts, expected intermediate block/repeat-table size, CPU/memory/disk
bound, software versions, restart behavior and no-paid-resource confirmation.
Small synthetic corruption controls may qualify parsing and restart behavior;
they are not a biological pilot and cannot narrow the full 526-taxon scope.

The extension is motivated by fungal genome-organization literature already
recorded in the [annotated bibliography](bibliography.md), especially the need
to distinguish transposable-element-rich compartments and structural variation
from protein-family turnover. It will be scheduled only after the current
structural-atlas and species-tree gates release their required inputs.

## Completion evidence

Completion requires public, source-bound gene-order/repeat/neighborhood tables;
reconciled anchor and synteny-block calls with uncertainty; branch/path mapping
across species-tree sensitivities; joint-model estimates with diagnostics and
multiple-testing correction; and figures/tables that preserve the relevant
missing/ambiguous dispositions. None of those results currently exists.
