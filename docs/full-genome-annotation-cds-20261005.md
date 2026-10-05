# Full assembly-to-annotation and CDS comparisons

This stage will compare original assembly DNA with all annotated coordinates
and available original CDS sequences across **526 entries**. It supports
annotation controls for sequence–structure coupling, duplication, domain
architecture and codon analyses. The full comparison has **not been launched
or queued**. Source freezing requires both full annotation and assembly-DNA
readbacks, with their actual original execution closures.

The software passed **19 offline checks** with actual original tool/native
exit zero. The checks are hand-calculated literal joins and complete small
source inventories, not a corpus pilot or biological qualification. A full
independent comparison reader remains to be prepared before interpreting
results. The overall project and all eight evolutionary aims remain unfinished.

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

Existing 519-taxon CDS-to-protein translation comparisons remain separate.
This stage compares DNA sequences and does not translate, infer selection,
or qualify a protein product as biologically correct.

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
will record exact qualified genome/target totals and resource estimates
before any expensive full launch. No GPU, prediction or paid provisioning is
part of this stage.

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

The blueprint is not a frozen full run plan. Its future source receipts are
unavailable at preparation. Run the freezer only after both full original
readbacks close; it binds every taxon/source and writes a distinct immutable
plan and estimates. There is no scheduler or queued full comparison.
