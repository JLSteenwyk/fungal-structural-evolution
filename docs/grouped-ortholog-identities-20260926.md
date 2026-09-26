# Full grouped ortholog identity audit

The completed small-family supplements checked small-family pairs but streamed
larger-family rows without validating their protein identities. The profile
native output contains 839,296,556 grouped rows. Before these tables can support
duplication–structure comparisons, every protein label must map to the stated
species and original protein family.

`scripts/audit_grouped_ortholog_identities.py` builds the complete source
protein/family lookup, then streams each native species table. It rejects
unknown species, same-species ortholog rows, unknown or wrongly assigned
proteins, malformed rows and repeated proteins within either side of a row.
Every native table must retain the SHA256 recorded by the completed supplement
scan. Source identity files are hashed before and after the audit. Original
outputs remain unchanged.

The audit counts rows and the product of the two grouped side lengths. These
are **directed incidences**, not unique ortholog pairs: cross-row duplication,
reciprocal equality, missing larger-family pairs and reconciliation semantics
are not established by this check. Biological orthology and duplication
interpretation remain separate requirements. No potentially billions-sized
flat pair table is materialized.

The existing native four-taxon, 18-gene software fixture passes: 20 grouped
rows and 22 directed incidences. Corrupted rows with an unknown family,
same-species target, missing protein, repeated protein or malformed columns
are rejected. The complete CLI also passed, and all output checksums and totals
were independently read back. Fixture evidence and its immutable execution
plan are `metadata/grouped_ortholog_identity_fixture_20260926.json`,
`metadata/grouped_ortholog_identity_fixture_plan_20260926.json` and
`metadata/grouped_ortholog_identity_fixture_completed_20260926.json`.
The first parser check predates the additional end-of-scan input hash check;
the CLI fixture and production plan pin the final script.

The full profile scan is now running under
`fungal-grouped-ortholog-identities-profile-20260926.service`. Launch identity
and configuration are recorded in
`metadata/grouped_ortholog_identity_profile_launch_20260926.json` and
`metadata/grouped_ortholog_identity_profile_plan_20260926.json`.

```bash
python scripts/audit_grouped_ortholog_identities.py --plan metadata/grouped_ortholog_identity_profile_plan_20260926.json
```

This command has been launched; do not run another copy. The output directory
is `results/orthology/grouped-ortholog-identities-profile-20260926-v1`.
`state.json` records completed tables after their hashes pass; only the final
`receipt.json` indicates completion of all tables. The process is capped at
one CPU and 32 GiB RAM, with no swap, low scheduling priority and a 50-GiB
initial free-disk gate. Planned output is below 1 GiB. The 1–48-hour planning
range is uncalibrated, not a completion forecast. This scan uses no GPU or paid
resources. The second guide and any full-scan numerical readback remain pending.
