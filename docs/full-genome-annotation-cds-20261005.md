# Full assembly-to-annotation and CDS comparisons

The full producer has compared original assembly DNA, annotation coordinates
and all available original CDS sequences across **526 entries**, retaining
**5,927,745 source products** and the unchanged **5,815,847 representatives**.
Original producer session **18749** exited zero; complete native, source-hash
and invocation-journal closure is recorded. The separate original reader
**56533** has reconstructed every entry and all **5,923,039 CDS records**,
with actual original API/native zero, full source-hash verification and
complete invocation-journal closure. The compact handoff binds **10,604
source/output files** and publishes all 526 taxon denominators.
[Producer receipt](../metadata/full_genome_annotation_cds_20261005_v1.json) and
[closed original execution](../metadata/full_genome_annotation_cds_transport_20261005_v1.json).

The producer reports **5,849,816 records with one exact unmodified genomic
CDS candidate**, **66,481 without an exact match**, **2,663 exact candidates
with retained annotation exceptions**, **3,988 missing/ambiguous product IDs**,
**90 records whose candidates all require review**, and **one exact candidate
with multiple target records**. These are target-record dispositions, not
gene-copy counts, translation qualification or evolutionary findings. These
counts agree with the completed independent every-record reconstruction.
[Full completed handoff](../metadata/full_genomic_cds_completed_20261005_v1.json),
[526-row taxon table](../metadata/full_genomic_cds_taxon_dispositions_20261005_v1.tsv)
and [separate table check](../metadata/full_genomic_cds_taxon_table_check_20261005_v1.json).

Software qualification retains **19 producer** and **18 independent-reader**
offline controls, without a corpus pilot or source repair. The full stages
use four CPUs, 64 GiB RAM, no swap, 56 GiB address space and one BLAS thread.
All eight evolutionary aims, genome contamination/haplotig checks, taxonomy
qualification and adequate ancestral uncertainty remain incomplete.

## Complete scope

| Component | Scope |
| --- | --- |
| Assembly/annotation entries | All 501 fungal entries and 25 outgroups |
| Source products | All 5,927,745 normalized products |
| Representatives | Unchanged baseline of 5,815,847 products |
| Annotation features | All 60,917,860 indexed feature rows |
| Provisional ORFs | Both complete tables; the out-of-bounds record retained |
| NCBI CDS targets | All 519 publisher files; 5,847,336 records and 7,954,834,617 bases in original receipts |
| External publisher CDS targets | Four deposited CDS FASTAs |
| External provisional targets | Both previous ORF-derived CDS files, explicitly not independent publisher evidence |
| Creolimax | All candidates/products retained; transcript sequence is not substituted for a CDS target |
| Complete CDS inventory | 5,923,039 records and 8,023,872,320 nucleotide bases |
| Complete genomic inventory | All 526 verified genomes and 26,645,610,508 bases |

Existing 519-taxon CDS-to-protein translation comparisons remain separate.
This stage compares DNA sequences and does not translate, infer selection,
or qualify a protein product as biologically correct.

The frozen plan binds **7,946 sources**, all entries, unchanged representatives
and all 60,917,860 original feature rows. CDS targets comprise 5,847,336 NCBI
records, 59,116 records from four publisher FASTAs and 16,587 previously
extracted provisional-ORF records. Those ORF targets retain their shared genome
dependency; the one out-of-bounds ORF remains a source/product review record.
Creolimax has no qualified original CDS target. Target-record and protein-product
counts are different inventories and do not establish gene-copy counts.
[Complete immutable source plan](../metadata/full_genome_annotation_cds_plan_20261005_v1.json)
and [source-derived resources](../metadata/full_genome_annotation_cds_resources_20261005_v1.json).

## Comparisons and retained dispositions

Every original annotation feature and provisional ORF receives a genomic
sequence-ID/coordinate disposition. All source proteins, alternatives,
selection decisions, uncertain mappings and source-unlinked CDS references
remain in scope. Genomic sequence identifiers match exactly; no new fuzzy
aliases are introduced. The provisional ORF tables already retain their
previously resolved genomic identifiers and original review statuses.

CDS rows group by the original NCBI feature ID, explicit published transcript
Parent, or GTF transcript candidate. GTF gene-level associations remain
separate candidate references. Repeated references to one feature are not
concatenated as duplicate exons. Missing/ambiguous feature groups receive
review dispositions.

Single-sequence joins use strand-aware coordinate order. A complete explicit
`part=X/Y` annotation provides order when present. Conflicting/incomplete
orders, unknown strands and unordered multi-sequence joins remain unresolved.
Overlapping intervals retain their repeated bases; phases and annotation
exceptions are recorded without trimming, recoding or repair. Documented
circular virtual intervals can be extracted; multipart circular joins without
explicit order require review. NCBI documents circular coordinates, overlapping
CDS parts and annotation exceptions, which is why this workflow retains these
cases explicitly.
[NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

Every available original CDS target retains its identifier, full header,
length and sequence digest. It is compared against every associated genomic
candidate, preserving all matching indices. Multiple exact loci, multiple
target records, gene-level matches and annotation exceptions retain review
status. A mismatch does not trigger frame search, phase removal or replacement
of the deposited target. One exact candidate is sequence agreement, not proof
of expression, gene-copy identity, annotation correctness or absence of
contamination/haplotigs.

Four per-taxon artifacts retain every coordinate, genomic candidate,
original CDS comparison and source-product disposition. Original genomes,
annotations, CDS files and the representative baseline remain unchanged.
Creolimax needs additional coding-boundary evidence; its original mRNA FASTA
is not an interchangeable coding-DNA target. ORF-derived targets are useful
consistency checks but share their original genome dependency.

## Software controls and resource design

The 19 checks cover plus/minus order, nonzero phases, overlaps, both-strand
circular intervals, explicit/conflicting/absent part order, bounds, mixed
strands, exceptions, multiple matching loci, multiple target records and
unrepaired mismatches. Complete literal accounting checks all four source
modes, alternatives/unresolved products, the provisional out-of-bounds row,
and the absence of a qualified Creolimax CDS target.

The resource design is four CPUs, 64 GiB RAM, no swap, 56 GiB per-process
address space, one BLAS thread and a 128 GiB output allowance. Whole per-taxon
genomes and CDS candidate/target caches motivate the RAM allowance. The
uncalibrated planning range is 2–72 hours with 2–32 GiB estimated output;
seven-day CPU/wall and 16 GiB per-file safety caps are not ETAs. The freezer
has recorded exact qualified genome/target totals and resource estimates
before any expensive full launch. No GPU, prediction or paid provisioning is
part of this stage.

## Independent reader

The separate reader imports no comparison-producer or join helpers. It uses
its own binary FASTA decoder, IUPAC byte-table reverse complement, feature
reference grouping and explicit-part ordering. It replays every feature
coordinate, genomic candidate, original target record and source-product
mapping/selection disposition. All matching loci, alternatives, phase and
exception fields and missing-target reasons must agree. Python, gzip, SQLite,
hashing libraries and original source inventories remain shared dependencies.

Eighteen offline controls cover all four source modes, hand-calculated IUPAC
complements, repeated gene/transcript references to one physical exon,
unavailable genomes and the unqualified Creolimax target. Ten deliberate
changes to coordinates, candidate/target sequence hashes, phase, part order,
exceptions, matching indices, target headers, alternative selection and an
unresolved product are rejected even after rebinding output hashes. These
literal checks are software evidence; the full corpus now has separate
independent replay and original execution closure.
The original single-call tool/native execution and complete wrapper journal
close. The completed full reader used four CPUs, 64 GiB RAM, no swap, 56 GiB
address space, one BLAS thread and an 8 GiB output allowance. Its 2–72 hour
planning range is uncalibrated; seven-day safety caps are not an ETA.

- [Separate reader](../scripts/readback_full_genome_annotation_cds_v1.py),
  [offline controls](../scripts/check_full_genomic_cds_readback_cases_v1.py),
  [actual result](../metadata/full_genomic_cds_readback_fixture_20261005_v1.json)
  and [original execution closure](../metadata/full_genomic_cds_readback_fixture_transport_20261005_v1.json).
- [Reader preparation blueprint](../metadata/full_genomic_cds_readback_blueprint_20261005_v1.json)
  and [full reader plan freezer](../scripts/prepare_full_genomic_cds_readback_plan_v1.py).

## Original launch and paired full readback plan

Both source and reader preparation commands have actual original tool zero
exits. The producer stage's exact wrapper/native launch, source plan/configuration
hashes and current cgroup limits are recorded separately from the immutable
prelaunch plan. The original 22:43 UTC launch checkpoint recorded input-hash checking
before output creation. Both full executions are now closed against their
actual original tool terminals and complete journals. Producer native wall
time was 21 minutes 8 seconds; the reader took 14 minutes 59 seconds. These
are observed runs, separate from the broad prelaunch resource estimate.
Mutable counters alone never establish full execution or sequence agreement.

- [Original full producer launch](../metadata/full_genome_annotation_cds_launch_20261005_v1.json)
  and [exact native runtime checkpoint](../metadata/full_genome_annotation_cds_checkpoint_20261005_goal_2244.json).
- [Full independent reader plan](../metadata/full_genomic_cds_readback_plan_20261005_v1.json)
  and [resources](../metadata/full_genomic_cds_readback_resources_20261005_v1.json).
- [Original source-plan preparation result](../metadata/full_genome_annotation_cds_preparation_original_tool_20261005_v1.json)
  and [paired reader-plan preparation result](../metadata/full_genomic_cds_readback_preparation_original_tool_20261005_v1.json).
- [Runtime observer](../scripts/record_full_genomic_cds_checkpoint_v1.py) checks
  the immutable plan/launch/current native identity; it does not repeat the
  producer's complete original-source scan or restart a missing handle.

## Reproducibility and launch gate

- [Prepared source/resource blueprint](../metadata/full_genome_annotation_cds_blueprint_20261005_v1.json).
- [Source freezer](../scripts/prepare_full_genome_annotation_cds_v1.py),
  [coordinate/join logic](../scripts/genomic_cds_join_v1.py) and
  [full comparison producer](../scripts/audit_full_genome_annotation_cds_v1.py).
- [Offline controls](../scripts/check_full_genome_annotation_cds_cases_v1.py),
  [actual result](../metadata/full_genome_annotation_cds_fixture_20261005_v1.json)
  and [original execution closure](../metadata/full_genome_annotation_cds_fixture_transport_20261005_v1.json).
- [Full annotation registry](full-annotation-coordinate-registry-20261005.md)
  and [full assembly-DNA acquisition/readback](full-assembly-dna-20261005.md).

The software blueprints retain their original preparation-time gates. The
separate frozen source plan above now binds both completed original source
readers and every taxon/source. The full comparison reader runs only after
actual original full producer API/native zero, full journal and source/output
closure. Immutable plans retain their prelaunch states; launch/execution
records separately establish what has actually run.

## Completed handoff and remaining biological controls

The compact merger compares every per-taxon producer/reader disposition and
publishes all 526 rows without a third full corpus scan. It accounts for
60,917,860 original annotation features and 16,588 provisional ORF rows:
60,934,448 coordinate rows total. All annotation features are within their
original sequence bounds; one original provisional ORF remains out of bounds.
There are 69 circular multipart candidates lacking qualified order and
21 candidates with incomplete/invalid explicit part order, all retained.

The merger's first invocation used a literal total 400 rows too high and
exited one before writing a table or handoff. The failed source and actual
tool result are retained. The corrected total is the sum of the two unchanged
source inventories; neither full native execution failed and no scientific
tolerance or disposition changed.
[Retained reporting failure](../metadata/full_genomic_cds_merger_reporting_failure_20261005_v1.json)
and [compact merger](../scripts/assemble_full_genomic_cds_completed_v1.py).

Large genomic comparison outputs remain outside Git at
`results/genome-annotation-cds-20261005-v1/`; independent replay is recorded
at `results/genome-annotation-cds-independent-readback-20261005-v1/`.
Reproduction requires the pinned original genomes, annotation databases and
CDS files, with fresh output namespaces. Scripts, plans, resource estimates,
closed execution logs, original tool payloads and all source hashes are
versioned; complete public release of the large data remains outstanding.

Next connect these dispositions to the existing CDS-to-protein translation
checks and preserve explicit protein/isoform/copy uncertainty in codon,
duplication and sequence–structure analyses. Resolve assembly redundancy,
contamination and taxonomic exceptions before biological copy admission.
Exact genomic DNA agreement alone does not qualify those conclusions.
