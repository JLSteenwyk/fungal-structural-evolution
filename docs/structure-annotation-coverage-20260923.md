# Annotation-stratified structural availability

This full-data analysis measures the frozen AlphaFold catalog's availability
within each taxon separately for representative proteins with any Pfam
 gathering-threshold hit and proteins without such a hit. It covers all
5,815,847 proteins, preserving 1,319,513 available model links. Counts and
fractions describe this source snapshot before confidence qualification;
local ESMFold additions are excluded.

The primary calculation streams all annotation protein/query links and joins
the original catalog TSV by exact taxon, protein and sequence identity. It
requires every catalog protein to be consumed once and every source protein
to be counted. A second calculation uses SQLite grouping and the independently
verified structure/family bridge to reconstruct all per-taxon annotation strata.
The two complete tables must agree, and every taxon's totals must match the
original catalog coverage table. Sequence mismatches, missing scope or differing
aggregates terminate the run. Full-CLI fixtures verified the proportions for
a three-protein example and rejection of a mismatched sequence.

Outputs include `taxon_annotation_coverage.tsv` and a receipt with pooled
counts, fractions and source bindings. Pooled and taxon-level contrasts should
be retained separately because taxon sampling can confound pooled comparisons.
Proteins sharing sequence or family and related taxa are not independent
replicates. No significance test, causal effect, domain-loss call or structural
novelty claim is made. Annotation status means any retained raw gathering-
threshold hit, not a validated functional assignment or structural domain.

Script: `scripts/summarize_structure_annotation_coverage.py`.
Plan: `metadata/structure_annotation_coverage_plan.json`.
Output: `results/structures/annotation-stratified-coverage-20260923-v1`.
Unit: `fungal-structure-annotation-coverage-20260923.service`.

The process was confirmed live and recorded in
`metadata/structure_annotation_coverage_launch.json`. Resources are one CPU,
8 GiB RAM, no swap and no GPU. Output planning is 0.01 GiB; runtime planning
is an uncalibrated 0.25–4 hours. The source hashes come from the previously
audited catalogs/databases and are checked before and after the run. Results
are pending; no coverage contrast is reported yet.
