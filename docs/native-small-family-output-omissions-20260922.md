# Small-family ortholog output omissions

Two isolated synthetic runs of the installed OrthoFinder expose incomplete
small-family ortholog tables despite successful exit and completed HOG/statistics
outputs. These are software-output limitations, not evidence of biological
absence. The production reconciliation is still running and has not been
accepted as a complete orthology result.

## Reproduced behavior

The fixture contains four taxa and seven families totaling 18 genes: one
four-gene family, three three-gene families spanning one/two/three taxa, two
two-gene families spanning one/two taxa, and a singleton. Only the four-gene
family requires a supplied gene tree. No real protein sequences are used.

With descending family sizes, 22 of 24 expected directed ortholog pairs are
written. The two-gene/two-taxon family is absent from the ortholog tables, even
though it is represented in the root HOG table. Inspection of the installed
`comparative_genomics/orthologues.py::AllOrthologues` shows that the two-gene
branch constructs its pair but never appends it to `all_orthologues`.

Putting the singleton immediately after the four-gene family reduces output
to 12 of 24 expected directed pairs. The same routine breaks at its first
singleton, skipping the subsequent three-gene families. The native `OGsAll`
reader preserves input group order; it does not sort families by size.
Our expanded partitions preserve source family identities and contain
non-singleton families after the first singleton.

These fixtures use the installed CLI with `--save-space --no-fix-files`, check
that staged inputs remain unchanged, and hash all output files. They record
missing pairs explicitly rather than reporting a passing complete output audit.
The observed counts apply to these fixtures; they do not validate biological
orthology, loss or duplication.

## Full input exposure

The assessment streamed both complete staged cluster files, checked each
family's full sorted membership against its crosswalk and verified the cluster
hash against the staged-file manifest. Both partitions contain all 5,815,847
proteins. A separate comparison reproduced all exposed family identifiers and
pair counts from the earlier reconciliation workload tables.

| Input partition | Two-gene/two-taxon families | Later cross-taxon three-gene families | Total exposed families | Expected missing directed pairs |
| --- | ---: | ---: | ---: | ---: |
| Profile guide | 31,056 | 11,567 | 42,623 | 124,056 |
| MAFFT guide | 30,890 | 11,447 | 42,337 | 123,144 |

The first singletons are OG0044000 and OG0046181, respectively. The totals are
anticipated omissions based on the tested native behavior, not yet measurements
of the running production ortholog tables. Directed counts include both
orientations of each pair and are not independent biological observations.

## Required handling

Preserve the native outputs and their provenance. Do not patch installed code
or the active controller while reconciliation is running, and do not relabel
families or rerun expensive large-family inference solely to recover this
small-family output. Full output validation must explicitly check these exposed
families and distinguish ortholog-table coverage from HOG/tree coverage.

A separately versioned small-family output supplement can reproduce the native
intended two-/three-gene rules without changing larger-family inference. Before
use, it must be tested against complete synthetic expectations and compared
against native production rows so missing pairs are filled without duplicating
existing rows. Any supplemented statistics must be labeled as recomputed; raw
native totals cannot silently be treated as corrected. The supplement implementation is now tested and queued as described below;
production supplementation and full native output readback remain outstanding.

## Reproduction and evidence

Each synthetic invocation uses two search threads and one analysis worker,
with a 120-second fixture-only timeout and small disposable output. The full
exposure scan uses one CPU and streams one family at a time; output is a few
megabytes. These checks neither run GPU prediction nor incur paid resources.

```bash
python scripts/check_native_small_family_outputs.py \
  --output results/orthology/native-small-family-output-fixture-v2
python scripts/check_native_small_family_outputs.py --early-singleton \
  --output results/orthology/native-small-family-early-singleton-fixture-v1
python scripts/assess_small_family_output_exposure.py \
  --inputs results/orthology/expanded-reconciliation-inputs-complete-v2 \
  --output results/orthology/small-family-output-exposure-v1
```

Commands require fresh output directories. Archived evidence:

- `metadata/native_small_family_sorted_observation.json`
- `metadata/native_small_family_early_singleton_observation.json`
- `metadata/native_small_family_exposure_receipt.json`
- `metadata/native_small_family_exposure_workload_readback.json`


## Tested supplement and queued production execution

`scripts/supplement_small_family_orthologs.py` enumerates the intended
cross-species directed pairs for every two-/three-gene family, independent of
family ordering. It scans every completed grouped ortholog file and compares
its small-family rows with that expectation. Unexpected pairs, duplicate pairs,
missing taxon tables and ambiguous labels fail validation. Large-family rows
are streamed without expanding their potentially enormous pair products; they
remain outside this check's validation scope.

Only absent pairs are written to `supplemental_directed_pairs.tsv`. Each record
includes the family, native gene IDs, species and protein labels and explicit
supplement provenance. Native files and statistics are never edited. The
receipt binds source mappings, the scanned native tables and the supplement.
The CLI requires a completed native execution, verifies the execution plan and
source manifest, and rejects modified bound inputs. It supports the tested
single-token accession naming mode; other normalization modes require review.

The full CLI produced two missing directed pairs for the size-ordered fixture
and twelve for the early-singleton fixture. Independent checks compared every
supplement record against the original fixture's missing-pair list and proved
that native plus supplemental small-family pairs equal the expectation without
overlap. Corruption tests reject duplicate/unexpected pairs, a missing species
table, changed source IDs and an incomplete execution. Reproduce with:

```bash
python scripts/check_small_family_supplement.py
```

The production controller `scripts/advance_small_family_supplement.py` is queued
under `fungal-small-family-supplement-20260922.service`, bound to the exact active
native controller PID/start time. It requires successful completion of both
guides before scanning either production result. The pinned plan is
`metadata/native_small_family_supplement_plan.json`; the launch record is
`metadata/native_small_family_supplement_launch.json`. Outputs will be under
`results/orthology/expanded-small-family-supplements-v1/`.

Resources are capped at one CPU, 8 GiB RAM and zero swap, with a 50-GiB free-disk
gate and 1-GiB planned output allowance. A 0.5–12 hour execution allowance is
uncalibrated planning, not an ETA; runtime depends on final native ortholog-table
volume. No GPU or paid resources are used. Full gene-tree/HOG/duplication
readback, final integrated ortholog statistics and biological interpretation
remain separate requirements.
