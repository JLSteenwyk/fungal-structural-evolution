# Taxon inputs for matched structural comparisons

Every one of the 2,786,912 selected records can now be joined to endpoint taxa,
species-kernel row indices, original family labels and shared-identity family
components. The exported table contains 65,287 unique nodes: 38,500 targets and
26,787 backgrounds. The original selected-record table remains authoritative
for policies, scenarios and matching orientation; it is checksum-bound in the
new receipt. This normalized representation preserves all records without
materializing millions of copies of node metadata.

The focal species occurs in the background pair in 265,580 selected records.
These are repeated records across settings, not 265,580 independent events.
The inventory records uses of each target and background node across the full
selection grid. Those counts describe reuse, not effective sample size or
analysis weights. Within any fitted setting, reuse must be recomputed after
structural qualification.

Endpoint labels retain their original graph order; they do not claim a residue
correspondence or a cross-comparison endpoint alignment. Both target endpoints
belong to the focal species. Background endpoints retain their individual taxa.
All taxa map into the fixed 526-tip order used by the five verified species
distance kernels. Original family labels are retained alongside the previously
verified combined gene/model/sequence identity components.

`scripts/prepare_selected_taxon_inputs.py` checks source artifacts and exports
`results/orthology/selected-taxon-inputs-20260927-v1/selected_node_taxa.tsv`.
`scripts/check_selected_taxon_inputs.py` independently reconstructs every table
cell using dataframe joins and grouped selection counts, and checks focal-species
overlap and family/guide agreement across all selected records. Its proof is
`metadata/selected_taxon_inputs_readback_20260927.json`.

This is model preparation. It does not fit an effect, choose a covariance
function, establish independent family blocks, or assume structural distances
are additive on the species tree. Those statistical choices and their
validation remain necessary before interpreting the observed contrasts.
