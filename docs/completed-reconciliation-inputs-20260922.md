# Completing expanded gene-tree inputs

The independently repaired OG0000017 tree has 55,981 distinct tips and SHA-256
`0b7af9780b1ffce2fc6c4d5abaa02ba37347431b94cba0ca382656beb1a46ffe`.
It was installed in the original tree collection, but remained absent from
both frozen expanded reconciliation input sets.

On September 22, a CPU job began creating new physical copies in
`results/orthology/expanded-reconciliation-inputs-complete-v2/`, preserving
`expanded-reconciliation-inputs-v1` unchanged. The new inputs include the
repaired tree in each guide partition, after checking its installation and
native-parser audit receipts, unique tips, membership checksum, family
crosswalk, and destination. All other copied files retain their previous
checksums. Expected tree totals are 95,484 for the profile partition and
95,594 for the MAFFT partition, each retaining all 5,815,847 proteins.

`scripts/complete_expanded_reconciliation_inputs.py` performs staging and then
runs `scripts/readback_expanded_reconciliation_inputs.py` against every copied
file and every tree. The independent native-parser check compares all tips
with the actual expanded clusters, checks branch lengths and exact file
universes, and verifies isolated native paths. Its completion status now
distinguishes a complete tree collection from a collection with explicit
pending families. Neither status claims completed reconciliation.

Configuration and process identity are recorded in
`metadata/expanded_reconciliation_inputs_complete_plan.json` and
`metadata/expanded_reconciliation_inputs_complete_launch.json`. The unit is
`fungal-complete-reconciliation-inputs-20260922.service`. Resources are two CPU
equivalents, 16 GiB memory, no swap, about 6 GiB additional data, and a 30 GiB
free-disk gate. The planning runtime allowance is 0.25–3 hours. Existing local
resources incur no new charges.

At launch this stage is **running**, with completion unproven. The authoritative
completion artifact will be the new directory's `readback.json` with status
`passed_complete_expanded_input_readback` and empty pending lists. Reconciliation
remains gated on a separate execution/resource plan and appropriate species-tree
inputs. No runnable native `Log.txt` or reconciliation outputs are created by
this preparation stage.
