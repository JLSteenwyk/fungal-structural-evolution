# Independent duplication-alignment readback

Both primary and additional-reference alignment runs now have a queued full
readback. The primary scope is 412,800 pair/order/mask dispositions; the
reference scope is 130,164. These jobs wait for the exact corresponding
producer process and require its matching successful receipt. They have not
yet processed production alignments.

`scripts/readback_duplication_alignments.py` reconstructs the full directed
grid from frozen pair tables, checks canonical pair hashes, exact checkpoint
paths, checkpoint-manifest hashes, model identities, commands, source-plan
hashes and input-manifest bindings. Every expected disposition must appear
exactly once. It verifies successful and unavailable/error/timeout outcomes
against their recorded inputs and native exit information. Failed outcomes
are checked for internal consistency, not rerun or independently adjudicated.

For every successful alignment,
`scripts/duplication_alignment_numeric_readback.py` independently reads saved
native output without using the producer parser, checks all alignment
characters and gap markers against hashed PDB sequences, and reconstructs
residue correspondences. It checks original-position mappings, aligned counts,
sequence identity, coverage and native text versus stored metrics. An
independent least-squares rigid fit checks RMSD to the native two-decimal
rounding tolerance (0.00501 Å); identity is checked to three-decimal rounding
(0.000501). The rigid-fit helper was previously used for cross-clan readback.
TM-scores are checked against native text and normalization labels, not
independently reoptimized. Matched pLDDT summaries use rounded PDB confidence.

Each readback uses one CPU, 16 GiB RAM, no swap and a bounded 256-input PDB
cache. BLAS thread counts are set to one. Estimated output allowance is 1 GiB
and the uncalibrated planning interval is 0.5–24 hours after its producer
finishes. No GPUs or paid resources are involved. Plans and exact process
identities are in:

- `metadata/duplication_alignment_readback_plan_20260926.json`
- `metadata/duplication_alignment_readback_launch_20260926.json`
- `metadata/duplication_reference_alignment_readback_plan_20260926.json`
- `metadata/duplication_reference_alignment_readback_launch_20260926.json`

Run either with `python scripts/readback_duplication_alignments.py --plan
<plan-path>`. Output directories are recorded in those plans. A successful
receipt will be named
`passed_full_duplication_alignment_mapping_rmsd_identity_readback` and include
all disposition counts, the number of numerically checked successes, maximum
RMSD rounding error and the hash of the numeric table.

Validation passed:

```bash
python scripts/check_duplication_alignment_numeric_readback.py
python scripts/check_duplication_alignment_full_readback.py
python scripts/check_duplication_primary_alignment_readback.py
```

The numeric fixture uses deformed native structures in both orders with full
and sparse confidence masks, checks nonzero reconstructed RMSD, and rejects
altered metrics and residue mappings. Completed-handoff fixtures exercise
both primary and two-source reference modes, successful native alignments
and excluded masks. Removing a disposition while updating its manifest hash
still fails the full-grid readback.

These checks do not establish structural novelty, ancestry, selection,
reference orthology, domain orientation reliability or a duplication effect.
Statistical analysis and interpretation remain separate requirements.
