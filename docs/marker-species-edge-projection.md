# Mapping marker branches to the species-tree guides

The completed projection maps all 18,397 paired-marker edges in the combined
ESMFold cohort and 26,758 in the expanded AlphaFold cohort to both provisional
526-tip guides. The cohorts are kept distinct, not treated as independent
biological replicates. This is topology correspondence, not a test of elevated
structural change.

For each full-guide edge, its bipartition is restricted to the marker's sampled
taxa. Empty sides are discarded and complementary splits are canonicalized.
Several full-guide edges can become the same pruned split: the marker edge then
covers a collapsed full-guide path, and its estimated change cannot be assigned
to a single constituent branch. A marker split absent from the pruned guide is
retained as discordant. Edge IDs encode canonical full-taxon splits and are
comparable between the two guides.

Of 45,155 marker edges, 20,399 map to the same single full edge in both guides;
9,943 map to the same collapsed path; 3,471 match both guides but their full-edge
mapping differs; 368 match only one guide; and 10,974 disagree with both.
Terminal edges inflate the apparent single-edge coverage. Internal-edge results
are therefore reported separately:

| Internal marker-edge outcome | Combined ESMFold | Expanded AlphaFold |
| --- | ---: | ---: |
| Same single full split in both guides | 1,975 | 4,600 |
| Same collapsed path in both guides | 1,245 | 1,753 |
| Matches both, mapping differs | 624 | 719 |
| Matches one guide only | 148 | 220 |
| Discordant with both | 5,073 | 5,901 |
| Total internal edges | 9,065 | 13,193 |

These guides are unrooted homogeneous-model trees without support estimates.
Even a unique agreement does not establish a reliable rooted clade assignment,
reconciliation, branch duration or acceleration. Gene-tree error, incomplete
sampling, discordance and model sensitivity remain relevant. The supported
mixture trees and remaining taxon/marker sensitivities must be incorporated
before final lineage tests. No marker rate or branch length has been divided
among the projected full edges, and no physical-displacement additivity is
assumed. The Gamma-optimization correction leaves these marker topologies
unchanged; this projection does not reuse its superseded numerical branch rates.

## Independent readback

An independent routine repeatedly pruned each full guide with Bio.Phylo for all
213 cohort-marker combinations (426 pruned trees). It reconstructed retained
split sets and checked all 90,310 projection rows, including compatibility and
path-size dispositions. For all 67,994 compatible rows, the pruned branch length
matches the sum of the identified original guide-edge lengths. This verifies
the topological projection and its path accounting, not the biological guide
or any evolutionary rate interpretation. The readback used one CPU core with
a 2-GiB planning allowance and a 0.01–1-hour runtime estimate.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/project_marker_edges_to_guides.py --guides results/phylogeny/profile-guide-audit-v1 results/phylogeny/mafft-guide-audit-v1 --fit-audits results/phylogeny/paired-fit-audit-esmfold-combined-v1 results/phylogeny/paired-fit-audit-gdm-expanded-v1 --fits results/phylogeny/paired-marker-fits-esmfold-combined-v1 results/phylogeny/paired-marker-fits-gdm-expanded-v1 --output results/phylogeny/marker-species-edge-projection-v1
OPENBLAS_NUM_THREADS=1 python scripts/audit_marker_edge_projection.py --projection results/phylogeny/marker-species-edge-projection-v1 --guides results/phylogeny/profile-guide-audit-v1 results/phylogeny/mafft-guide-audit-v1 --output metadata/marker_species_projection_readback.json
```

Use fresh outputs. Full edge identities, projection rows, marker summaries and
guide sensitivity tables remain in results/phylogeny/marker-species-edge-projection-v1.
Receipts, readback, resource plan and the terminal/internal cohort summary are
versioned under metadata/marker_species_projection_*.

Next steps are to repeat this mapping against supported species-tree variants,
retain path/discordance uncertainty in branch-level models, and assess
sequence–structure change within homologous proteins and resolved regions.
The current count of uniquely mappable edges is an eligibility diagnostic,
not a completed acceleration or clade-association analysis.

## Numerical branch frame

`results/phylogeny/projected-gamma-branch-frame-v1` attaches four paired G4 branch estimates (AA, 3Di AF exchangeabilities, 3Di AF with empirical frequencies, and 3Di LLM exchangeabilities) to every one of the 45,155 marker edges. The 852 fitted trees contribute 180,620 branch lengths. It uses the revised ESMFold Gamma baseline, including the reviewed optimization correction, and the audited expanded AlphaFold Gamma fits. Alignment site counts and fit-warning counts remain explicit.

Every row retains both full-guide edge sets and the prior mapping disposition. A full-edge ID is assigned only when both guides identify the same single edge. Collapsed paths and discordant edges keep their marker-level estimates without allocating change among constituent species branches. The two prediction-source cohorts remain separate.

Unique internal-edge coverage is 188 edges for ESMFold (median eight markers, 71 edges with at least ten) and 139 for AlphaFold (median twenty markers, 81 edges with at least ten). These are coverage diagnostics, not thresholds for statistical significance or independent replication. Terminal edges are summarized separately in `metadata/projected_gamma_branch_coverage_summary.tsv`.

The independent readback reconstructs all tree splits through undirected graph edge deletion, checks all 180,620 lengths, preserves all 90,310 guide projections and checks each unique-edge coverage group. This validates assembly of the numerical frame. Branch lengths are expected substitutions in the fitted alphabet per aligned site, not physical displacement or change per unit time. The table is not an acceleration test; supported-guide, FreeRate, family, coverage and phylogenetic-dependence analyses remain required. The procedure is implemented in `scripts/assemble_projected_branch_frame.py` and `scripts/audit_projected_branch_frame.py`; the source hashes and readback are versioned under `metadata/projected_gamma_branch_frame_*`.
