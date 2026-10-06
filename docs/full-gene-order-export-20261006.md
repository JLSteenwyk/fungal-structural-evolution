# Full immutable gene-order export

This is the first executable full-panel input for genome-context analyses. It
will export one gzip table per frozen manifest taxon from the completed
annotation-coordinate registry. Native annotation `gene` features retain their
source scaffold, coordinates, strand, feature identity, within-scaffold order,
neighbor gaps and linked source/selected products. The two taxa represented by
verified ORF coordinates remain explicit provisional rows and are never
converted into gene anchors.

The [resource plan](../metadata/full_gene_order_resources_20261006_v1.json)
uses one read-only SQLite input at a time: two CPUs, 8 GiB memory, zero swap and
a 32-GiB output allowance. Its 1–48-hour range is intentionally uncalibrated.
It is based on the full 526-taxon receipt census, not a pilot subsample.

`build_full_gene_order_v1.py` rejects a changed manifest/registry taxon set,
receipt disagreement, missing database, duplicate registry receipt, pre-existing
output, and any taxon whose emitted coordinate rows fail the source receipt's
expected gene or ORF count. `check_full_gene_order_v1.py` separately validates
every compressed table hash, schema, taxon identity, coordinate ordering,
within-scaffold ordinal and bidirectional neighbor-gap arithmetic.

The output is a coordinate/context table only. Product links are native source
mappings, not reconciled orthology or synteny anchors. No repeat, collinear
block, rearrangement, copy-number or structural-evolution result follows from
this export.
