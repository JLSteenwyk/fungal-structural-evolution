# Full coding and taxonomic genetic-code context

The full source audit compares every inherited CDS translation test code to
the pinned taxonomic nuclear and mitochondrial assignments. It covers all
**501 fungal entries plus 25 outgroups**, **5,927,745 original protein products**,
**5,815,847 representatives** and **5,923,039 original CDS targets**. All
111,898 alternative products and 67,810 separate derived-evidence products
remain in scope. This prepares genetic-code review for family codon analyses;
it does not replace any original translation, CDS target or protein sequence.

## Source interpretation

The [NCBI primary schema](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt)
documents separate nuclear and mitochondrial code identifiers and inheritance
flags in the common first 13 `nodes.dmp` fields. The current schema text is
cached with its URL/time/hash; it does not update the historical taxonomic
snapshot. The five additional new-taxdump fields, including plastid and
hydrogenosome assignments, are not interpreted by this audit.

The [NCBI genetic-code reference](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
describes CDS-level translation-table choices. The original audits' missing
explicit-code policy remains an explicitly labeled table-1 assumption.
Taxonomic assignment or exact translation under a test code does not establish
the actual code or the genomic compartment of a particular source sequence.
Missing mitochondrial assignment (`0`) remains unspecified, not code 1.

The 161,820,772-byte frozen `new_taxdump.tar.gz` and the already verified
526-entry identity table supply the source context. Every selected species
path is rechecked against its original taxonomic path. Nuclear and mitochondrial
inheritance are traced independently to their explicit ancestor, retaining
raw values and inheritance flags. Assembly catalogue identifiers supply an
additional context for 521 entries; five publisher depositions retain no
invented NCBI assembly context. Earlier hybrid, uncertain-name, species-complex
and strain cautions remain attached. This is not accepted species delimitation
or a current-taxonomy claim.

## Frozen assignments and record classifications

The prepared context retains nuclear assignments for **516 code-1 entries,
eight code-12 entries, one code-6 entry and one code-26 entry**. All inherited
raw values agree with their resolved assignments. Mitochondrial assignments
are separate: 474 code-4 entries, 15 code-3, two code-16, 24 code-1, five code-5
and six unspecified. Neither distribution establishes sequence compartment.

Every original target retains its ordinal, CDS/protein identity, DNA hash,
original joined status, inherited translation status, test code and code
provenance. Relationships distinguish nuclear agreement, agreement with both
compartments, mitochondrial-only agreement requiring compartment verification,
disagreement with both assignments, and missing/unspecified test codes.
No code is substituted or chosen by whether it improves translation.

Every original product retains its sequence hash, representative decision,
all original target contexts, model availability and separate derived code
tests. The four publisher projections and Creolimax extraction preserve their
original code-1 tests, original outcomes and shared genome dependency. They
do not become original publisher CDSs. Alternatives retain the same baseline
assignment scope as the structural union.

## Completed execution and results

Twenty literal inheritance/identity/classification controls pass with actual
original software zero and complete wrapper/source closure. The first reader
fixture failed because its parent output directory had not been created. Its
V1 source, actual original/native exit 1, whole invocation journal and partial
synthetic outputs are retained. V2 adds only that fixture-directory creation;
the producer helper, independent reader and acceptance conditions are unchanged.
V2 passes the complete synthetic archive/identity/target/product integration
and three semantic corruption controls. These controls invoke semantic checks
directly; no hash-rebinding exercise or real taxon pilot is performed.

Source preparation completed under actual original session 58837/native zero
with whole source/wrapper/journal closure. The plan freezes **1,103 bindings
and 2,119,231,430 input bytes** before full launch. The producer completed all
526 entries under actual original session 3496/native zero, with **2,165
source/output bindings** and invocation start/end closure. Its measured runtime
was 2 minutes 13 seconds. The separate reader launched only after this closure;
actual original session 98890 completed all 526 entries with native/API zero.
Its whole wrapper messages, source/output bindings and invocation start/end
close. Its measured runtime was 2 minutes 32 seconds. It independently parses
the raw archive, resolves every ancestor and reconstructs
every original target/product record without producer/helper imports. Python
libraries and the original coding/taxon evidence remain shared dependencies.

Both paths verify **20,478 products** with a test code different from both
snapshot assignments. This includes **17,211 representatives** with an available
model and original genomic/inherited strict-translation agreement:

| Frozen entry | All products with code disagreement | Original exact agreement plus a model |
| --- | ---: | ---: |
| Candida africana, F241526 | 2,850 | 427 |
| Ascoidea rubescens, F54195 | 6,787 | 6,326 |
| Candida viswanathii, F5486 | 10,841 | 10,458 |

Their snapshot nuclear assignment is code 12. These are source discrepancies
requiring genetic-code, sequence-compartment and taxonomy review, not proved
annotation errors. Another **235 original-exact modeled representatives** match
only the mitochondrial assignment and require compartment verification. The
completed handoff merges **2,179 source/output bindings** and publishes all
526 taxon rows. Every public table identity, inherited-code ancestor and
status/model count passes a further check against the independent full replay.
The original 2,952,711 modeled representatives with genomic/strict-translation
agreement remain unchanged; 2,935,265 match the nuclear assignment (including
90,886 whose nuclear and mitochondrial assignments coincide). That agreement
does not resolve compartment, biological genetic code or selection eligibility.

Both full stages use four CPUs, 64 GiB RAM, no swap, a 56 GiB per-process AS
limit and one BLAS thread. The 0.1–12 hour/0.5–16 GiB planning ranges are
uncalibrated; seven-day CPU/wall safety bounds are not an ETA. Preparation used
one CPU/16 GiB RAM/12 GiB AS. No GPU prediction, translation, optimization,
paid service or biological inference is launched.

## Reproduction and evidence

The [full plan](../metadata/full_genetic_code_context_plan_20261006_v1.json),
[complete pinned taxonomic contexts](../metadata/full_taxonomic_genetic_code_context_20261006_v1.json),
[original producer launch](../metadata/full_genetic_code_context_launch_20261006_v1.json)
and [original reader launch](../metadata/full_genetic_code_context_readback_launch_20261006_v1.json)
record the sources and actual executions. The
[completed handoff](../metadata/full_genetic_code_context_completed_20261006_v1.json),
[all 526 taxon rows](../metadata/full_genetic_code_context_taxon_dispositions_20261006_v1.tsv)
and [complete table check](../metadata/full_genetic_code_context_taxon_table_readback_20261006_v2.json)
record the final verified counts. The
[source-preparation closure](../metadata/full_genetic_code_context_preparation_transport_20261006_v1.json)
and [retained fixture failure](../metadata/genetic_code_context_reader_software_failure_20261006_v1.json)
keep success and failure distinct. Large per-record gzip files and synthetic
fixtures remain outside Git history.

The following record the original commands. Preserve existing receipts and use
fresh source-bound output/plan/receipt namespaces for reproduction; the source
preparation script contains fixed original source paths. Recreating a historical
taxonomy download requires the recorded snapshot hash, not a silent new download.

```bash
python scripts/prepare_full_genetic_code_context_v1.py \
  --plan metadata/full_genetic_code_context_plan_20261006_v1.json \
  --context-output metadata/full_taxonomic_genetic_code_context_20261006_v1.json
python scripts/build_full_genetic_code_context_v1.py \
  --plan metadata/full_genetic_code_context_plan_20261006_v1.json \
  --receipt metadata/full_genetic_code_context_20261006_v1.json
```

Genetic-code adjudication, compartment/annotation review, full independent
DNA translation, family codon alignment/divergence, copy/homology checks,
confidence/predictor calibration and accepted phylogenetic uncertainty remain
required. All eight evolutionary aims are incomplete; GPU prediction stays paused.
