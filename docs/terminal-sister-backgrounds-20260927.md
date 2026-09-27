# Terminal sister-pair background inventory

The duplication analysis needs comparison groups before testing an association
between duplication and structural divergence. This inventory begins that work
across all 70,307 profile-guide and 70,412 MAFFT-guide resolved gene trees.
It is running; no completed background or duplication-effect test is claimed.

For every internal node with exactly two terminal children, retain both genes,
their taxa, family, node and parent identifiers, reported duplication flags for
the node and parent, both terminal sequence branch lengths and their sum.
Join each gene to the frozen structure bridge and retain four availability
states: two distinct models, identical model, one model, or no models. Missing
models do not remove pairs from the inventory. Internal or multifurcating
sister groups are outside this terminal-pair design; this is not an inventory
of all possible ortholog pairs or all evolutionary events.

Each pair receives one descriptive label, with reported duplication taking
precedence: reported duplication, same-taxon unreported, or cross-taxon
unreported candidate. **An unreported duplication is not proof of speciation.**
Cross-taxon candidates can retain deeper duplication ancestry and losses.
Parent duplication flags are retained separately rather than silently filtering
them away. The two guide outputs are dependent sensitivity alternatives.

Next steps are full independent reconstruction of the inventory, cross-guide
gene-pair comparison, and ortholog-membership checks before matching candidate
backgrounds to duplicated pairs by family and sequence divergence. Domain
architecture, length, confidence, taxon sampling and shared ancestry also need
to enter eligibility, matching or downstream models. Cross-species comparisons
and same-species duplicates differ in evolutionary context; matching alone will
not justify a causal duplication interpretation. No structural outcome is used
to select this initial pool.

Run `scripts/inventory_terminal_sister_backgrounds.py --plan
metadata/terminal_sister_background_plan_20260927.json`. Inputs and imported
helpers are checksum-pinned. Output is under
`results/orthology/terminal-sister-background-inventory-20260927-v1`; each guide
has a TSV and the final receipt records exact counts and hashes. The output
directory must be new. An interrupted run must be diagnosed before any new
version is launched; no partial TSV is a completed inventory.

Known-tree fixtures in `scripts/check_terminal_sister_background_cases.py`
passed canonical gene orientation, path lengths, reported node/parent flags,
all model availability states, same-taxon and multifurcation distinctions,
multiple cherries, and rejection of invalid branches, duplicated labels and
unknown taxa. These fixtures do not replace the full production readback.

The launch reserves one CPU, 8 GiB RAM, no swap, and a 2 GiB output allowance;
the planning range is 0.1–4 hours, not a calibrated ETA. Free disk was about
12 TiB and available RAM about 942 GiB before launch. The service enforces CPU
and memory limits. No GPU predictions or paid resources are used. Process
identity and the source-bound plan are recorded in
`metadata/terminal_sister_background_launch_20260927.json`.
