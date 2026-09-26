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

## Matching MAFFT-guide audit launched

The matching MAFFT scan has now launched using the same fixture-tested script
and resource caps. Its 526 completed native tables total 41,159,858,967 bytes;
the source snapshot records 839,250,464 grouped rows. The launch preflight
verified all table paths, the unchanged checker/parser and approximately
12.28 TB of free disk. Each guide has independent inputs and outputs, so the
two read-only scans run concurrently, capped at two CPU equivalents and 64 GiB
RAM in total. At launch the profile scan was live and using about 634 MiB RAM.

The MAFFT plan and live process identity are
`metadata/grouped_ortholog_identity_mafft_plan_20260926.json` and
`metadata/grouped_ortholog_identity_mafft_launch_20260926.json`.
The service is `fungal-grouped-ortholog-identities-mafft-20260926.service`,
writing `results/orthology/grouped-ortholog-identities-mafft-20260926-v1`.

```bash
python scripts/audit_grouped_ortholog_identities.py --plan metadata/grouped_ortholog_identity_mafft_plan_20260926.json
```

This command is already running; do not duplicate it. Both complete receipts,
their output readbacks and comparison of the guide-specific results remain
outstanding. The identity audit still does not validate ortholog-pair semantics
or justify interpreting native directed incidences as independent observations.

## Reproducible validation command

The previously recorded manual fixture checks now have a standalone runner:

```bash
python scripts/check_grouped_ortholog_identities.py --output results/orthology/grouped-ortholog-identity-reproducibility-20260926-v1
```

This run passed. It rechecks the native fixture source hashes, rejects six
malformed/incorrect row cases, executes the full audit CLI, independently
enumerates the small fixture's grouped pairs to verify each table count, and
confirms rejection of a changed provenance pin and an existing output path.
It uses the existing four-taxon/18-gene software fixture, runs on one CPU in
under a second, and does not rerun native phylogenetic inference. All runner
artifacts were subsequently rehashed; the receipt is archived in
`metadata/grouped_ortholog_identity_reproducibility_20260926.json`.
Use a fresh output path to reproduce it. The two live production scans and
their pinned checker were unchanged.

## Large-field parser recovery

The first profile scan terminated with exit code 1 after 164 complete tables
(264,544,191 rows). Python's default CSV field limit rejected a 174,200-byte
gene-list field in F181123.tsv, line 40,897, family OG0000017. The file still
matches the completed native snapshot. Reproducing the failing parse and then
allowing a 16-MiB field verified all identities in the offending row and its
12,443 directed incidences. This is a parser-size issue, not an established
biological expansion or an identity mismatch.

The original script, failed service and partial output are preserved.
`scripts/audit_grouped_ortholog_identities_large_fields.py` sets the CSV field
limit and calls the unchanged producer. Its full fixture run passed with
unchanged 4-taxon/18-gene/20-row/22-incidence totals; output hashes were read
back. The actual-row and fixture checks are archived in
`metadata/grouped_ortholog_large_field_check_20260926.json` and
`metadata/grouped_ortholog_large_field_fixture_20260926.json`.

The profile scan was relaunched only after confirming terminal failure, using
fresh output `results/orthology/grouped-ortholog-identities-profile-20260926-v2`
and the original resource caps. Its plan and launch identity are
`metadata/grouped_ortholog_identity_profile_large_fields_plan_20260926.json`
and `metadata/grouped_ortholog_identity_profile_large_fields_launch_20260926.json`.
The plan pins both the wrapper and original producer. This restarts the audit
scan, not native reconciliation. Earlier state counts are not reused as
completed output artifacts. The service is
`fungal-grouped-ortholog-identities-profile-large-fields-20260926.service`.
The original MAFFT audit was still active at this recovery checkpoint and was
not altered. Full audit completion remains outstanding.

The original MAFFT audit subsequently terminated with the same CSV limit
failure. Its offending OG0000017 row independently passed the larger-field
identity check (12,443 directed incidences), recorded in
`metadata/grouped_ortholog_large_field_mafft_check_20260926.json`. After terminal
failure was confirmed, it too was relaunched through the tested wrapper with
fresh `grouped-ortholog-identities-mafft-20260926-v2` output and unchanged resource
caps. Its plan/launch records are
`metadata/grouped_ortholog_identity_mafft_large_fields_plan_20260926.json` and
`metadata/grouped_ortholog_identity_mafft_large_fields_launch_20260926.json`;
service `fungal-grouped-ortholog-identities-mafft-large-fields-20260926.service`.
Both v1 scans remain failed attempts; neither partial count establishes full
validation. Both v2 scans are running and their final receipts remain pending.

## Count-table readback prepared

`scripts/readback_grouped_ortholog_counts.py` requires a completed identity
audit and independently sums both its family and taxon tables. It checks
unique keys, integer counts, the complete native taxon-table set, recorded
native hashes against the pinned snapshot, source protein count, and the
plan/producer/source/count-file hashes. It does not rescan the native ortholog
rows or independently validate their per-family assignments or pair semantics.

The four-taxon fixture passed (18 proteins, 20 rows, 22 directed incidences);
the result is archived in
`metadata/grouped_ortholog_counts_fixture_readback_20260926.json`. Temporary
copies with refreshed artifact hashes still rejected unequal family/table
totals, a duplicate family key, and a missing taxon table. The unchanged
temporary copy passed. These checks did not alter production data.

After each v2 scan writes its final receipt, run the matching command below.
Neither production readback has been executed yet:

```bash
python scripts/readback_grouped_ortholog_counts.py --plan metadata/grouped_ortholog_identity_profile_large_fields_plan_20260926.json --output metadata/grouped_ortholog_counts_profile_readback_20260926.json
python scripts/readback_grouped_ortholog_counts.py --plan metadata/grouped_ortholog_identity_mafft_large_fields_plan_20260926.json --output metadata/grouped_ortholog_counts_mafft_readback_20260926.json
```

Each readback uses one CPU, streams source hashes and protein counts, and
loads only the two small count tables. Planned resources are under 2 GiB RAM,
under 1 MiB output, and minutes of runtime (uncalibrated); it does not repeat
the approximately 41-GB native-table scan. Full identity and biological
orthology validation remain distinct from this numerical consistency check.

## Profile-guide scan completed

The corrected profile-guide service exited successfully (exit code 0), and
its final receipt is archived in
`metadata/grouped_ortholog_identity_profile_completed_receipt_20260926.json`.
All 526 native tables passed the identity scan: 839,296,556 grouped rows,
1,123,279,884 directed row incidences and 62,742 families with rows, using
the full 5,815,847-protein source partition. These incidence counts are not
deduplicated ortholog-pair counts.

The production count-table readback command above also passed; evidence is
`metadata/grouped_ortholog_counts_profile_readback_20260926.json`.
Both count-table totals agree with the receipt, all 526 recorded native hashes
match the pinned snapshot, and all checked source/producer/artifact bindings
remain intact. The MAFFT-guide scan is still running at this checkpoint.
Reciprocity, cross-row duplicates, missing larger-family pairs and biological
orthology remain unvalidated, as do their use in duplication–structure tests.

## MAFFT-guide scan completed

The corrected MAFFT-guide service subsequently exited with code 0. Its
archived receipt is
`metadata/grouped_ortholog_identity_mafft_completed_receipt_20260926.json`:
526 tables, 839,250,464 grouped rows, 1,123,035,370 directed row incidences,
62,903 families with rows and 5,815,847 source proteins. Its production
count-table and provenance readback also passed, archived in
`metadata/grouped_ortholog_counts_mafft_readback_20260926.json`.
Both full identity scans and their count readbacks are now complete; the
earlier running states above are historical. The pair-semantics and biological
limitations described above still apply to both guides.
