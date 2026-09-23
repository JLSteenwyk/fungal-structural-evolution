# Hierarchical orthogroup output verification

The profile-guide native reconciliation has completed, and exact resolved-tree
membership checks passed for 70,307 trees. Hierarchical orthogroup (HOG)
outputs require separate validation before duplication and domain analyses.
The full profile output contains 525 node tables totaling 6,381,538,832 bytes.

`scripts/audit_hog_output_identities.py` is running against every table under
`fungal-profile-hog-identities-20260922-v2.service`. The plan
`metadata/profile_hog_identity_audit_recovery_plan.json` pins all tables, native-stage
receipt, source species/protein/family maps, labeled species tree and audit
implementation (532 files). It checks:

- The complete expected set of internal-node tables and source-taxon headers.
- Valid, unique and consecutive node-specific HOG identifiers.
- Every emitted protein's exact taxon and original family identity.
- Containment of every assigned species within the named species-tree clade.
- No repeated gene assignment within one hierarchical level, while permitting
  expected reuse across different levels.
- Stable output files during reading and unchanged source hashes afterward.

The verifier independently reads source membership and does not call the
native HOG inference code. The installed native writer uses all source-species
columns at every level, so containment is tested on nonempty cells rather than
by restricting table headers. Native two-/three-gene HOG output uses a separate
routine that continues past singleton families; the previously observed
small-family orthologue-pair omissions must not be assumed to apply to HOGs.

The output will be `results/orthology/profile-hog-identities-20260922-v2`, with
per-node row/assignment counts and hashes. An assignment count summed across
levels is not a number of distinct proteins. The checks do not establish
complete assignment of all eligible genes, correctness of ancestral group
membership, gene-tree parent-clade labels, hierarchical nesting, rooting,
orthology or duplication events. These remain further requirements, as does
verification of the still-running MAFFT-guide reconciliation.

A hand-constructed three-taxon fixture permits a gene to recur across levels
and rejects wrong-family assignments, genes outside the node's species clade,
within-level duplicate genes, unknown protein identities and a missing table.
Run it with `python scripts/check_hog_output_identities.py`. The fixture passed;
full biological-output verification is still running.

Pre-launch resources are one CPU, 32 GiB RAM, no swap, 0.01 GiB output, and an
uncalibrated 0.1–4 hour planning range on existing local resources. No GPU or
paid service is used. The native output and running reconciliation are read-only.

The first audit attempt terminated with exit status 1 at Python's default
131,072-byte CSV field limit. Large same-species gene-list cells exceed that
limit. Its output directory and launch plan are preserved; no successful
scientific verification was claimed. The recovery raises the parser limit
to a bounded 16 MiB and uses a fresh output directory and service. An added
160-kilobyte cell fixture verifies that semantic duplicate checks still run
past the old parser limit. The original implementation remains in Git at
commit `2451379`; recovery plans and launch records identify the corrected run.

## Completed full identity audit and root disposition

After the September 23 reboot, the read-only audit completed in fresh output
`results/orthology/profile-hog-identities-20260923-v3`. All 525 node tables
passed: 9,277,304 HOG rows and 93,034,008 gene assignments summed across levels.
The latter is not a distinct-gene count. The completed receipt and all per-node
counts/hashes are archived as `metadata/profile_hog_identity_completed_readback.json`
and `metadata/profile_hog_identity_node_summary.tsv`.

A separate full-data source/root comparison now accounts for missing root
assignments using `scripts/resolve_root_hog_gene_disposition.py`:

| Disposition | Genes |
|---|---:|
| Source proteins | 5,815,847 |
| Source singleton families (outside this root comparison) | 500,113 |
| Source proteins in nonsingleton families | 5,315,734 |
| Assigned to 183,956 root HOGs | 5,274,605 |
| Missing from root HOGs | 41,129 |
| Missing and in native phylogenetically misplaced lists | 35,187 |
| Missing without that native flag | 5,942 |

Every native flagged gene is absent from the root HOG table. However, the
flagged-gene lists do **not** fully explain the missing set. All 5,942 remaining
genes belong to families of four or more genes; those families are covered by
the previously passed full resolved-tree membership readback. Their absence
from root HOGs therefore does not imply absence from the resolved gene trees.
The native HOG writer explicitly skips leaves carrying its `X` flag. Its other
ancestral-group and duplication rules still need examination before assigning
a cause to the remaining omissions. No root-completeness, biological gene-loss,
contamination or horizontal-transfer claim is justified by these counts.

The new analysis checks every emitted root protein against the original
species/protein/family maps, rejects repeated or wrong-family assignments,
retains all unassigned nonsingleton genes and compares their exact identities
with all 519 native flagged-gene files. All missing genes span 5,688 families
and 521 taxa. Per-gene details remain outside Git at
`results/orthology/profile-root-hog-disposition-20260923-v1/missing_genes.tsv`;
its checksum, counts and per-taxon summary are archived in
`metadata/profile_root_hog_disposition_receipt.json` and
`metadata/profile_root_hog_taxon_disposition.tsv`.

The plan pins 525 input/implementation files, uses one CPU and 32 GiB RAM
without swap, allows 0.1 GiB output and planned 1–30 minutes. It completed in
approximately 31 CPU seconds. Four full-CLI fixture cases distinguish an
exactly explained omission from an unexplained omission and reject unknown
flagged identities and wrong source-family assignments. Reproduce with
`python scripts/check_root_hog_disposition.py`; full-data reproduction requires
a fresh output path in `metadata/profile_root_hog_disposition_plan.json`.

## Resolved-tree placement of every unflagged root omission

The September 23 parent-clade trace completed for all **5,942 unflagged missing
root genes**, spanning **1,712 families and 486 taxa**. All lie outside the
emitted root HOG parent clades. In these affected families, every one of the
**24,116 emitted root HOG memberships** equals the descendants of its named
resolved-tree parent node after removing native flagged genes. Thus the
unflagged missing genes are not dropped from within those represented clades.
Every omitted gene's immediate parent also has at least one emitted HOG below
it; none requires climbing an additional ancestral edge to reach such a join.

These are properties of the native output partition, not independent evidence
that the ancestral grouping or duplication assignments are biologically correct.
The native rules that selected these parent clades still need semantic review.
Keep the omitted genes in the source inventory and distinguish them from
biological absences in subsequent duplication, loss and structure comparisons.
No missing gene was reassigned and no native output was changed.

Reproduce with `scripts/trace_root_hog_omissions.py --plan
metadata/profile_root_hog_omission_trace_plan.json` using a fresh output path.
The pinned plan uses one CPU, 16 GiB RAM, no swap and no GPU. Full results remain
at `results/orthology/profile-root-hog-omission-trace-20260923-v1`; the archived
receipt is `metadata/profile_root_hog_omission_trace_receipt.json`.
`python scripts/check_root_hog_omission_trace.py` passed exact placement and
multi-edge join fixtures and rejects wrong membership, nested root clades,
missing targets and already assigned targets. A separate full-table readback
matched every output identity to the source missing-gene inventory, checked
uniqueness and counts, and summarized join distances; its record is
`metadata/profile_root_hog_omission_trace_readback.json`. That readback checks
identities and table consistency, not an independent reconstruction of trees.

## Native writer mechanism reproduced

A minimal, read-only reproduction now exercises the installed `HogWriter`
directly, using species tree `(0,1)N0` and gene tree
`(0_a,((0_b,1_b)n2,(0_c,1_c)n3)n1)n0`. With `n1` assigned a duplication
at `N0`, the writer emits root HOGs at `n2` and `n3`, omitting `0_a`.
With that duplication flag removed, the same topology yields one root HOG
at `n0`, including all five genes. No misplaced-gene flags are used.

The mechanism combines the writer's suppression of ancestral HOGs above a
duplication with its immediate return for leaf nodes: the isolated leaf does
not receive its own fallback HOG. The reproduction assigns the reconciliation
features explicitly; it does not validate the upstream inference of those
features. It demonstrates how an unflagged gene can be omitted without being
absent from the source or tree. The production placement results are consistent
with this mechanism, but attributing all 5,942 cases to it still requires
reconstructing and checking their actual native reconciliation state.

Run `.cache/envs/orthofinder/bin/python scripts/check_native_hog_leaf_omission.py`.
The archived result `metadata/native_hog_leaf_omission_reproduction.json`
records the implementation hash, both outcomes and scope. It uses negligible
CPU/memory, opens no native output files and changes neither the installed
package nor production group assignments. Future loss analyses must retain
explicit unassigned states rather than converting these missing assignments
to inferred losses.

## Production-family replay

A deterministic selection of the smallest affected family chose **OG0010309**
(seven genes). Reclassifying its existing resolved tree with the installed
native reconciliation code, then running the native HOG writer, reproduced
both emitted root parent clades (`n2`, `n5`) and their exact memberships.
The same gene, `O400682_XP_011409098.1`, was omitted; no genes were marked as
misplaced in this replay. The reconstructed duplication at gene-tree node
`n1` maps to species-tree node `N0`, matching the mechanism above. Two other
duplications (`n3`, `n4`) map to native species 30. The native root-duplication
support is 0.0038022813688212928; this is the software's species-overlap
fraction, not a bootstrap probability or independent evidence of an ancient
duplication.

Reproduce with `.cache/envs/orthofinder/bin/python
scripts/replay_root_hog_omission_case.py`. The result is archived at
`metadata/native_root_hog_production_case_replay.json`, with source and code
checksums. This small CPU-only diagnostic completed in seconds and did not
modify native outputs. It starts from the already resolved production tree:
initial rooting and tree resolution are not independently rerun. It confirms
one real instance rather than establishing the cause of all omissions or the
biological correctness of the duplication assignment. The full affected-family
replay and guide sensitivity remain outstanding.

## All affected families: native replay running

The full **1,712-family** replay launched on September 23 under
`fungal-profile-root-hog-native-replay-20260923.service`, with original native
species and gene identifiers. It covers every family containing an unflagged
root omission, without a family-size exclusion. For each resolved tree it reruns
native classification and the HOG writer, compares exact parent-clade
memberships, misplaced-gene identities and missing identities, and records
the immediate parent state of every target gene. Mismatches and per-family
errors are preserved; a passing receipt requires all families to match.

The pinned plan `metadata/profile_root_hog_native_replay_plan.json` binds 57
source and implementation files, including the installed native Python sources.
`metadata/profile_root_hog_native_replay_launch.json` records the verified live
process. Limits are one CPU, 32 GiB memory, no swap, no GPU or paid resources;
output planning is 1 GiB and runtime planning is an uncalibrated 1–48 hours.
Available memory exceeded 800 GiB and disk space exceeded 12 TiB at planning.
Small families run first; every completed family's result is flushed to
`results/orthology/profile-root-hog-native-replay-20260923-v1/family_replays.jsonl`.
The run uses a fresh directory and leaves production outputs untouched.

This is a replay from resolved trees, not independent initial rooting or tree
resolution. Biological duplication/loss interpretation and alternative-guide
sensitivity remain separate requirements, even if all memberships match.

## Final omission accounting queued

The full native replay remains active. At the 1,668-family checkpoint all
families matched, and all 5,432 target genes processed so far had immediate
parents mapped to `N0` with an `N0` duplication at or below that parent. This
is a partial observation, not a result for the remaining families.

A final accounting controller now waits for that exact native-replay process.
It requires a passing full receipt and pinned output hash, checks the complete
source family/gene grid and flagged-gene counts, and matches every target's
parent label to the independent resolved-tree placement trace. It also requires
the same number of affected root HOG memberships as the earlier full placement
analysis. It emits an explicit algorithmic disposition for each gene; any other
parent state is retained for review rather than forced into the proposed rule.

Script: `scripts/summarize_root_hog_replay.py`.
Plan: `metadata/profile_root_hog_replay_disposition_plan.json`.
Output: `results/orthology/profile-root-hog-replay-dispositions-20260923-v1`.
Unit: `fungal-root-hog-replay-dispositions-20260923.service`.
The controller is confirmed live and recorded in metadata. Limits are one CPU,
2 GiB memory and no swap; planning allows 0.01 GiB output and 1–5 minutes after
replay completion. Full-CLI fixtures passed exact identity/count/parent-state
accounting and rejected a conflicting independently recorded parent label.

This check accounts for the native algorithm's output on the existing resolved
trees. It does not validate the ancestral rooting, duplication timing or
biological homology, nor interpret unassigned genes as losses. Overall
reconciliation interpretation and the other-guide comparison remain unfinished.
