# Full assembly-to-annotation and CDS comparisons

This stage will compare original assembly DNA with all annotated coordinates
and available original CDS sequences across **526 entries**. It supports
annotation controls for sequence–structure coupling, duplication, domain
architecture and codon analyses. The complete source plan is now frozen after both full annotation and
assembly-DNA readers passed with actual original execution closures. The full
comparison launched under original session **18749** with four CPUs, 64 GiB RAM
and no swap. At 22:43 UTC it is verified live while checking source hashes;
no comparison rows have been completed yet. The paired full reader plan is
prepared, not launched or queued before actual original producer closure.

The software passed **19 offline checks** with actual original tool/native
exit zero. The checks are hand-calculated literal joins and complete small
source inventories, not a corpus pilot or biological qualification. A separate full
comparison reader is prepared and passes **18 offline controls**. Its full
source plan is prepared after completed genomic source verification; full
reader execution remains gated on the original full comparison producer. The overall project and all eight evolutionary aims remain unfinished.

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

Four per-taxon artifacts will retain every coordinate, genomic candidate,
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
literal checks are software evidence; the full corpus has not been compared.
The original single-call tool/native execution and complete wrapper journal
close. The prepared full reader uses four CPUs, 64 GiB RAM, no swap, 56 GiB
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
prelaunch plan. At 22:43 UTC no per-taxon outputs exist: input-hash checking
precedes creation of the output root. Mutable counters never establish full
execution or sequence agreement.

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
