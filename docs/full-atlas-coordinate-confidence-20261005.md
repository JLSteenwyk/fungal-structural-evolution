# Full source coordinate and residue-confidence audit — October 5, 2026

The whole atlas now enters coordinate-content and residue-confidence auditing.
The **complete 2,961,055-model source universe** is running under eight CPU
workers, not a marker subset or separate pilot. The inputs preserve both
predictors, original settings and the unlinked annotated alternative product.
They are the independently closed
[full availability union](full-prediction-atlas-union-20261005.md).
No additional structures are predicted and GPU prediction remains paused.

| Source | Models attempted in the full plan | Inventory residue positions | Maximum length |
| --- | ---: | ---: | ---: |
| AFDB | 2,935,733 | 1,152,706,482 | 2,699 |
| ESMFold | 25,322 | 12,221,520 | 1,024 |
| Combined | 2,961,055 | 1,164,928,002 | 2,699 |

These are source-model residue positions, including predictor overlaps; they
are not independent biological observations or accepted high-confidence
residues. The full source inventory was scanned to verify all three dimensions
before launching. Representative availability stays 3,019,669/5,815,847 (51.9%),
with 2,796,178 proteins lacking either model. All 526 taxon denominators remain
in the upstream dataset, including missing proteins and source-specific gaps.

For every source model the audit uses the unchanged qualified
`extract_domain_coordinates.load_atoms` validator. It rehashes the original
coordinate file and verifies the canonical polymer sequence/hash/length,
atom/residue identity, one chain/model, unique atoms, finite coordinates and
occupancy/confidence bounds, and complete C-alpha coverage. Unsupported sources
are retained with explicit rejection dispositions. Changed immutable source
bytes cause a stage failure, not a confidence exclusion. Missing backbone
positions are recorded explicitly and never removed from the protein inventory.
This strict supported-representation policy is not proof that a rejected protein
cannot fold or has no function.

Accepted representations export exact float64 C-alpha coordinates and confidence
values, an exact sequence byte array and a backbone-missing mask. Compressed NPZ
members are stored in 1,481 tar shards, with at most 2,000 source models per
shard, avoiding millions of new filesystem entries. No extra coordinate or score
rounding is introduced. All atom coordinates remain in the bound original files.
Every member, array and source has a checksum, with full source-specific counts.

Inclusive pLDDT cutoffs 50, 70 and 90 are descriptive sensitivity masks. They
are not whole-protein admission criteria or calibrated accuracy estimates.
Original predictor/configuration identity is retained; scores are not used to
choose a winner between predictors. Coordinate-field and catalog confidence
means are both reported, with their difference; representations can differ in
score precision, especially converted local PDB files. PAE is explicitly marked
unevaluated. Full native structural-alphabet context confidence, PAE availability,
directional matrix qualification and predictor calibration remain separate work.

Two bounded literal controls have actual original zero exits:

- The producer control verifies sub-PDB-rounding coordinates, sequence,
  scores 49.99/70, missing backbone, rejection of nonfinite coordinates,
  out-of-bounds confidence, residue mismatch, incomplete C-alpha and multiple
  chains, source-bound tar serialization, and hard failure for changed bytes.
- The independent-reader control verifies literal coordinates and exact
  50/90 score boundaries, scalar threshold recounts, complete job/archive
  membership and rejection of wrong counts, modified self-consistent arrays
  and changed immutable source bytes.

These are software/representation checks with literal fixtures, not a corpus
pilot or source accuracy calibration. The native producer source is unchanged
after its qualified control. A separate reader reconstructs arrays from original
CIF columns without importing producer measurement code, verifies every full
source/job row and archive member, and rederives rejection eligibility. It uses
the same CIF lexical parser and canonical residue map, which remain shared
dependencies. Exact recorded rejection causation is not independently adjudicated.

Both controls finished in their first tool call, with no tool session created.
The first read-only evidence checker incorrectly required an optional manager
resource-completion message, absent for the subsecond run, and failed before
writing a result. Its source, actual failure payload and audit are preserved.
A new V2 checker verifies the actual original tool exit and complete original
wrapper header/native terminal, records the manager completion count as zero,
and invents neither a tool session nor a manager record. The successful native
fixtures were not repeated. See the
[producer control](../metadata/full_atlas_coordinate_profile_fixture_20261005_v1.json),
[its exact execution evidence](../metadata/full_atlas_coordinate_profile_fixture_transport_20261005_v2.json),
[reader control](../metadata/full_atlas_profile_readback_fixture_20261005_v1.json),
[its exact execution evidence](../metadata/full_atlas_profile_readback_fixture_transport_20261005_v1.json)
and [preserved observational failure](../metadata/full_atlas_coordinate_profile_fixture_transport_failure_audit_20261005_v1.json).

The original full audit, invocation `3406863b6435439a92f2abb25d21947d`,
completed with actual exit zero in 40,963.44 seconds. Its immutable execution
closure is
[`metadata/full_atlas_coordinate_profiles_execution_20261005_v1.json`](../metadata/full_atlas_coordinate_profiles_execution_20261005_v1.json).
The original wrapper start and terminal records, plus the recorded external
tool-session identity, are independently bound in
[`metadata/full_atlas_coordinate_profiles_transport_20261006_v1.json`](../metadata/full_atlas_coordinate_profiles_transport_20261006_v1.json).
The external wait terminal itself was not retained; the closure uses the
original bounded execution receipt and matching systemd terminal instead of
inventing one.
The
[immutable plan](../metadata/full_atlas_coordinate_profiles_plan_20261005_v1.json),
[bounded resources](../metadata/full_atlas_coordinate_profiles_resources_20261005_v1.json),
[exact launch](../metadata/full_atlas_coordinate_profiles_launch_20261005_v1.json)
and [original initial tool payload](../metadata/full_atlas_coordinate_profiles_original_initial_tool_20261005_v1.json)
preserve the source and execution identities.

Before launch, resource estimates specify eight CPUs, 64 GiB memory, no swap,
one BLAS thread and 56 GiB per-process address-space limits. Input coordinates
previously verified total approximately 1.105 TB. Arrays occupy 39.61 GB before
compression, NPZ/tar overhead and metadata. The 750 GiB output allowance,
2 TiB initial free-space requirement and 1 TiB emergency reserve fit the more
than 10 TiB available storage. No paid resource or new download is used.
Scaling the prior full four-worker all-atom extraction gives approximately
14.6 hours at eight workers; parsing, compression, source-model length and
I/O differences are uncalibrated. The planning range is 8–168 hours; the
168-hour hard cap is not a finish ETA.

The completed producer emitted all 1,481 planned shards. It reports 2,935,733
valid AFDB models and 25,322 valid ESMFold models, with no rejected model or
missing-backbone residue. This covers 1,164,928,002 source-model residues and
does not qualify the atlas for structural biology, provide a prediction-accuracy
estimate, or replace independent readback. Its final receipt binds the saved
output and exact totals; independent reconstruction of source atoms, arrays,
archive membership and dispositions remains pending.

The independent reader is prepared and its literal control passes. It has
**not been launched**; its full immutable input plan must bind the completed
producer and actual original transport first. It will check every original
source/job record, all original atom identities, exact exported sequence and
coordinate/confidence arrays, dtypes/shapes/digests, scalar threshold totals,
archive members and rejection eligibility. This requires another full source
read and separately estimated resources before launch.

The recorded execution configuration gives the complete bounded launch command.
Its native stage is:

```bash
python scripts/audit_full_atlas_coordinates_v1.py \
  --plan metadata/full_atlas_coordinate_profiles_plan_20261005_v1.json \
  --receipt metadata/full_atlas_coordinate_profiles_20261005_v1.json
```

The data and exact pinned inputs must be available; existing outputs are refused.
Use the recorded wrapper/cgroup limits for a fresh reproduction. Coordinates,
jobs, tar profiles and per-model dispositions remain outside Git history.

The separate full [domain extraction](completed-domain-extraction-20261005.md)
and its atom/archive readback have also completed. Full biological fits remain
zero. Accepted phylogenetic/reconciliation/dating frameworks, calibrated sequence–structure
effects, adequate ancestral uncertainty, all eight evolutionary aims and
publication deliverables remain required and unfinished.
