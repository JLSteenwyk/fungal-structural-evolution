# Exploratory structural similarity groups

The first clustering stage includes every model in the two frozen marker-model
snapshots: 13,153 AlphaFold and 5,121 ESMFold models, 18,274 total. Each receives
a unique alias linked to its original source, model ID/version, sequence hash,
coordinate path/hash, length and mean CA pLDDT. Identical sequences from different
prediction sources remain separate model records. Taxon and marker mappings
remain available in the pinned original snapshots; model clusters do not remove
species-specific variation from evolutionary inputs.

Foldseek runs its 3Di+AA alignment mode with sensitivity 7.5, up to 2,000
prefilter hits per query, E-value cutoff 1e-5, 80% aligned coverage in both
query and target and an approximate alignment-normalized TM-score threshold of
0.5. Greedy length-based single-step clustering and reassignment are requested.
The full command, installed executable hash/version and source receipts are
pinned in `metadata/frozen_structure_clustering_config.json`; option semantics
were checked against the installed executable's `easy-cluster -h` output.

pLDDT below 70 is masked for seeding only. This does not remove every uncertain
residue from alignment or enforce the six-residue/PAE gates used by paired
phylogenetic analyses. Clusters are exploratory similarity groups, not
orthogroups, validated homolog groups, functions or novel folds. Prefilter hit
limits may reduce recall. Representatives and members require explicit edge
validation; cluster membership does not imply all member pairs meet thresholds.
Threshold, prediction-source, domain-boundary and confidence sensitivities remain
required before interpretation. Full-chain coverage can miss valid domain-level
relationships, which require their own atlas analysis.

```bash
python scripts/cluster_frozen_marker_structures.py --output results/structural_clusters/frozen-marker-models-v1
```

The run is already active. It uses eight CPU threads, a 32 GB prefilter split
memory limit, 48 GB total memory planning and 100 GB output planning on the
existing host. The initial 1–48 hour forecast is broad pending production
timings; no paid infrastructure was provisioned. Inputs are symlinks to verified
source files, and temporary databases are retained for validation. The wrapper
rejects an existing output directory; inspect any partial failure before recovery
rather than starting a duplicate or overwriting results.

All 18,274 aliases and symlink targets passed readback, and the producer verified
every coordinate checksum. Foldseek database construction has started; there is
no completed clustering result yet. This is the full frozen **marker-model**
set, not the project's entire fungal proteome atlas.


## Initial clustering completed; membership verified

The run completed with 2,249 similarity groups across all 18,274 model records:
1,276 singleton groups and 302 groups containing both prediction sources.
Every Foldseek input database alias occurs exactly once. Independent readback
verified the complete output partition, representative self-membership, unique
annotated members, every retained provenance field, and agreement between raw
and annotated assignments. Source-specific taxon/marker links remain untouched.

These are algorithmic group counts, not numbers of orthogroups, novel folds or
independent evolutionary events. Alignment criteria have not yet been
independently checked for each representative/member pair. In particular,
confidence seeding does not certify all aligned residues. Edge validation,
threshold sensitivity, confidence/source analysis and domain-level clustering
remain pending. Large retained intermediate databases allow those follow-ups.
