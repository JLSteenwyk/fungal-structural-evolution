# OG0000017 source-annotation and repeat-context audit

## Purpose

OG0000017 contains 55,981 genes across the sampled taxa. One focal assembly,
*Austropuccinia psidii* GCA_000469055.2 (taxon `F181123`), supplies 46,787
members. This audit retains every one of the assembly's 133,973 selected
proteins and measures annotation context before the family is used in a
duplication, domain-turnover, or structural-evolution result.

## Inputs and method

The audit joins native orthology identifiers to the original annotation-coordinate
registry, assembly metrics, exact UniProt sequence-match table, and all 64
completed Pfam-chunk files. It records protein and coding length, gene span,
number of CDS segments, scaffold length and gene count, adjacent-gene gaps,
hypothetical-product labels, exact UniProt matches, Pfam presence, and a
descriptive flag for repeat-associated terms in Pfam names/descriptions.

Repeated amino-acid sequences are measured instead of being discarded. The
complete run has 133,953 distinct source sequence hashes among 133,973 products;
39 products occur in 19 non-unique sequence-hash groups.

## Result and interpretation boundary

The successful v3 receipt reports 46,787 focal-family proteins and 87,186 other
focal-assembly proteins. The family cohort has 19,764 proteins with the
descriptive repeat-keyword Pfam flag, compared with 6,046 in the other cohort.
All focal products have a `hypothetical protein` CDS label in this source
annotation. These observations identify a material annotation/repeat-risk
sensitivity for this family. They do not classify any record as a transposable
element, change an orthogroup membership, establish duplication, or support a
structural or ecological evolutionary conclusion.

The v1 and v2 executions are retained as failed, bounded software records. V1
incorrectly rejected duplicate source sequence hashes; v2 reached a strict TSV
writer schema error. V3 corrects only those software issues and supplies the
accepted complete result.

## Reproduction and validation

Run `scripts/audit_og0000017_annotation_repeat_context_v1.py` with the exact
paths and checksums recorded in
`metadata/orthology_OG0000017_annotation_repeat_context_20261006_v3.json`.
Then run `scripts/check_og0000017_annotation_repeat_context_v1.py` against the
receipt and `protein_context.tsv`. The checker verifies the receipt hash, table
schema, complete 133,973-product census, unique native/protein IDs, and the
46,787/87,186 cohort split.
