# Genome-architecture resource census

This receipt-level census supports planning the full 526-taxon genome-context
extension. It reads only the completed coordinate-registry and assembly-statistics
receipts; it does not reopen the 64.5-GB registry, call repeats or synteny, or
launch a biological analysis.

The immutable [census receipt](../metadata/genome_architecture_resource_census_20261006_v1.json)
records one disposition for every frozen manifest taxon. It checks that each
assembly-statistics receipt, where available, is verified and names exactly the
same assembly accession as the manifest. A deliberate altered-accession control
was rejected before output was written.

## Measured planning inputs

All 526 taxa have completed coordinate-registry receipts: 60,917,860 annotation
feature rows, 23,456,867 CDS rows, 5,983,922 gene rows and 5,815,847 retained
representative proteins. The registry itself occupies 64,512,704,512 bytes.

Matched, verified assembly-statistics receipts exist for 519 taxa, covering
26,327,586,003 assembly bases and 551,852 scaffolds. The remaining seven taxa
are explicitly marked `not_available`; their missing statistics must remain a
stratum or be recovered before any analysis that needs assembly-contiguity
covariates. The receipt preserves all mapping modes, including the two
provisional ORF-coordinate taxa.

For immutable flat-text planning only, the census uses conservative stated
assumptions of 512 bytes per gene-order row and 384 bytes per context row. They
bound a gene-order table at 2,977,713,664 bytes, a per-gene context table at
2,233,285,248 bytes, and an immediate left/right neighbor table at
4,466,570,496 bytes (11,631,694 directed neighbor rows at most). These are
storage-planning bounds, not observations of actual output sizes.

## What remains before a full launch

The census does not bound repeat-annotation or collinear-block intermediates,
runtime, peak memory, or disk required by their software. Those bounds need
the selected versioned repeat and synteny implementations, reconciled anchors,
and a restartable execution design. The full launch remains gated on those
plans as well as the active structural-atlas and species-tree sensitivity
workflows. No genome-context, rearrangement, or structure-association result
is reported here.

Reproduce with `scripts/census_genome_architecture_resources_v1.py`, using the
paths and hashes recorded in the receipt.
