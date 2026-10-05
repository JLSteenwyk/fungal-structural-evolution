# Full annotation coordinates and protein-product registry

This stage indexes the original annotations and all protein products across
**526 entries: 501 fungi and 25 outgroups**. It supports genome-to-CDS checks
and controls for annotation artifacts in duplication, domain architecture and
sequence–structure analyses. The registry does not establish gene correctness,
assembly quality, expression, gene-copy counts or evolutionary events.

The full producer and independent reader have completed across all **526
entries, 60,917,860 features, 5,927,745 source products** and the unchanged
5,815,847 representatives. All 111,898 alternatives and 16,588 provisional ORF
rows remain available. The databases total **64,512,704,512 bytes**.

Original producer session **9848** and reader session **29353** returned zero.
Every source feature/nonfeature line, Parent/candidate reference, original
mapping/selection row, protein sequence digest and representative record was
replayed. Complete native/manager journals and all declared bindings closed;
the compact handoff merges **3,721 bindings**. This verifies the index and its
source custody, not annotation correctness or assembly-to-CDS agreement.
[Completed registry handoff](../metadata/full_annotation_coordinate_completed_20261005_v2.json),
[full independent replay](../metadata/full_annotation_coordinate_readback_20261005_v1.json)
and [original reader execution closure](../metadata/full_annotation_coordinate_readback_transport_20261005_v1.json).

Producer completion was at 21:21 UTC; reader completion was at 21:37 UTC and
reader execution closure at 21:42 UTC on October 5. Earlier 21:08/21:34 partial
checkpoints remain unchanged. The reader's prepared plan still records its
prelaunch state; actual launch/completion are separate immutable evidence.

## Complete scope and unchanged baseline

| Source | Scope |
| --- | --- |
| NCBI GFF3 | All 519 exact assembly/annotation versions |
| Published external annotations | Four transcript GFFs and the Creolimax GTF |
| Sanchytrid coordinates | Both complete provisional ORF tables |
| Normalized source proteins | All 5,927,745 products |
| Selected representatives | The unchanged 5,815,847-product baseline |
| Alternative products | All 111,898 excluded products retained in the registry |

The NCBI receipts enumerate **60,434,823 features**, including **23,306,261
CDS segments**, in 1,183,724,940 compressed annotation bytes. The full plan
binds 2,643 original source, script and control files; per-taxon source hashes
also bind normalized proteins, original gene maps, representative FASTAs and
selection decisions. Original manifest rows and annotation provenance remain
in the plan. No taxon, product or uncertain mapping is removed to obtain a
cleaner result. Unique-species sampling and biological source qualification
remain separate unresolved project requirements.

Each taxon has a SQLite database with:

- Every original feature line, source line number, genomic sequence identifier,
  feature type, coordinates, strand, phase, score and decoded attributes.
- Original nonfeature lines, directives and any embedded FASTA tail.
- Every feature-to-Parent reference, including CDS rows without an ID.
- Explicit CDS product references or direct transcript/GTF candidate links,
  retaining references outside the source proteome and products without a
  direct CDS candidate.
- Every source protein identifier, header, length, sequence digest, original
  gene-mapping row, selection decision and representative membership.
- Separate provisional ORF rows with their original status, including the
  known Sanchytrium coordinate-out-of-bounds record.

Multipart feature IDs remain repeated rows. Attribute lists split on literal
commas before percent decoding; encoded commas remain inside their values.
Circular/virtual coordinates, overlaps, phase, partial flags, `part`,
`transl_except`, `exception`, and genetic-code attributes are retained. Bare
published GTF labels are explicitly marked as unparsed labels and preserved
in the original line. No phase trimming, coordinate repair, product reassignment
or promotion of ORFs to gene models occurs.

NCBI documents multipart IDs, origin-spanning coordinates, overlapping CDS
segments and annotation exceptions; these cases require review rather than
automatic trimming or coordinate rejection.
[NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

## Candidate associations and independent checks

NCBI CDS references use explicit `protein_id`; external GFF references use
direct transcript Parent candidates. GTF transcript links and gene-level
candidate links remain distinct. Gene-level candidates require additional
review; indirect transcript links are not silently invented. The entire
feature/Parent graph is available for subsequent resolution. The registry
does not concatenate, translate or compare CDS sequences with genomic DNA.

The independent reader imports no producer decoder. It reconstructs every
original feature and nonfeature line with separate percent decoding, checks
every Parent/candidate reference, and replays all source/mapping/decision and
representative records with a separate FASTA parser. Counts, fields, sequence
digests, exceptions and source hashes must match. Python, SQLite and Biopython
remain shared dependencies; this is not a wholly independent software stack.

Both versions of the offline controls passed **16 checks** with actual
original tool and native exit zero. They cover all four source modes,
multipart/phase/circular cases, escaped commas and translation exceptions,
alternatives/unresolved/outside-source products, and deliberate changes to
coordinates, phase, feature IDs, Parent links, product references, selection,
gene mappings, missing features and exception attributes. Corruption checks
rebind the altered database hash so source hashing alone cannot produce a pass.
These are literal software cases, not a corpus pilot or biological benchmark.

V1 and its successful evidence remain unchanged. Prelaunch review identified
a possible race when counting transient SQLite journal files during another
worker's commit. V2 counts stable database paths, with journal/index scratch
included in the resource allowance. No full V1 corpus run was launched and no
annotation rule or numerical tolerance was changed.

## Resources and continuation gates

Estimates were recorded before full launch: **four CPU workers, 32 GiB RAM,
no swap, 28 GiB address space per native process and one BLAS thread**.
The output estimate is 20–128 GiB; the declared allowance is 256 GiB including
journal/index scratch. A periodic stable-database guard at 192 GiB leaves
scratch/in-flight headroom; this is not an atomic filesystem quota. A 16 GiB
per-file cap and 1 TiB emergency free-disk reserve also apply. The uncalibrated
planning range is 1–48 hours; seven-day CPU/wall caps are safety limits, not
ETAs. No GPU, new structure prediction, paid provisioning or charge is used.

The full independent reader requires the actual original producer API/native
zero, complete original journal and all source/output bindings. Those producer
gates closed before its full launch with four CPUs and 32 GiB/no swap. The [subsequent assembly-to-CDS stage](full-genome-annotation-cds-20261005.md)
also requires independently qualified full assembly DNA. Existing
519-taxon CDS-to-protein comparisons are retained; they do not replace genomic
DNA reconstruction. Annotation exceptions, circular joins, alternative
products and uncertain matches must remain explicit in that stage.

The full structural atlas, accepted phylogenetic frameworks, adequate
ancestral uncertainty, all eight evolutionary aims and final release remain
unfinished.

## Reproducibility

- [Full plan](../metadata/full_annotation_coordinate_plan_20261005_v2.json),
  [resources](../metadata/full_annotation_coordinate_resources_20261005_v2.json),
  [original launch](../metadata/full_annotation_coordinate_launch_20261005_v2.json)
  and [runtime checkpoint](../metadata/full_annotation_coordinate_checkpoint_20261005_goal_2109.json).
- [Source freezer](../scripts/prepare_full_annotation_coordinate_registry_20261005_v2.py)
  and [full producer](../scripts/build_full_annotation_coordinate_registry_v2.py).
- [Full reader](../scripts/readback_full_annotation_coordinate_registry_v1.py),
  [reader plan](../metadata/full_annotation_coordinate_readback_plan_20261005_v1.json)
  and [reader resources](../metadata/full_annotation_coordinate_readback_resources_20261005_v1.json).
- [Literal controls](../scripts/check_full_annotation_coordinate_registry_cases_v2.py),
  [result](../metadata/full_annotation_coordinate_fixture_20261005_v2.json)
  and [original execution evidence](../metadata/full_annotation_coordinate_fixture_transport_20261005_v2.json).

Large databases remain outside Git under
`results/annotation-coordinate-registry-20261005-v2/`.
