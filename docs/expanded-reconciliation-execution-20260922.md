# Full expanded reconciliation execution

The execution controller launched September 22 with both complete expanded
partitions: 658,183 profile-guide families and 658,522 MAFFT-guide families,
each containing all 5,815,847 proteins across 526 taxa. The existing rooted
homogeneous guides define conditional sensitivities. Supported-mixture
species-tree alternatives, rooting uncertainty and biological interpretation
remain separate requirements.

`scripts/run_expanded_reconciliation.py` first checks the complete independent
input readback and pinned native software. It makes new physical copies of
every input file, then invokes the installed OrthoFinder CLI with
`--from-trees`, the matching named species tree, `--no-fix-files`, and
`--save-space`. The profile run precedes the MAFFT run; each uses four analysis
workers and eight search threads. Existing validated input snapshots remain
outside the execution directories.

The exact four-worker command configuration passed a synthetic software
fixture: 13 genes, three families, all 42 directed ortholog pairs, one expected
terminal duplication, and unchanged input files. The missing-ID-tree fixture
again returned exit code zero without completed outputs. Consequently the
controller requires final HOG, resolved-tree, duplication, species-tree and
statistics files, correct species/gene totals, and unchanged copied inputs
before calling a stage execution-complete. Full output readback remains
mandatory afterward; execution completion alone is insufficient.

The installed native restart path hardcodes `save_space=False` for its first
inference call. With file fixing disabled, `--save-space` therefore suppresses
the later split into species-pair tables but leaves plain per-species TSVs.
The fixture confirmed this behavior. No compression savings are assumed, and
no installed library was modified. The source implementations are pinned.

The service `fungal-expanded-reconciliation-20260922.service` caps CPU use at
eight equivalents and memory at 128 GiB, disables swap, and allows 65,536 file
descriptors. Each stage requires 192 GiB available memory and 3 TiB free disk
before launch. A watchdog checks disk every 20 seconds and terminates only
the native process group if free space falls below 1 TiB. Its state then
requires review; it does not automatically restart incomplete inference.
No GPU or paid resource is used.

The output planning range is 100–2,560 GiB per guide, plus about 6 GiB of
physical input copies across both guides. The 0.5–21-day runtime range per
guide is an uncalibrated planning allowance, not an ETA. The measured workload
contains about 5.2 billion candidate cross-species protein pairs per guide;
these are not inferred ortholog counts or a forecast of grouped output size.
All families, including the annotation-sensitive large OG0000017 family, are
retained. Its concentration in one taxon must not be interpreted as a validated
biological expansion.

Plan: `metadata/expanded_reconciliation_execution_plan.json`.
Live state and isolated execution directories:
`results/orthology/expanded-reconciliation-execution-v1/`.
Fixture evidence:
`metadata/orthology_native_save_space_four_worker_fixture_{receipt,readback}.json`.
The first stage begins with input copying; native inference starts only after
those copies pass hash checks.
