# Alternative-copy tree review

This review investigates alternative gene copies for observations whose coding
segments overlap publisher FCS EXCLUDE/FIX/TRIM regions. These are exploratory
gene trees, not accepted ortholog replacements. The 33-marker queue includes
38 coverage-qualified alternative copies alongside the original selected
proteins and available baseline taxa. Trees use the same profile-alignment
coordinates, restricted amino-acid model selection and 1,000 SH-aLRT replicates.

## First completed marker: 4764044at2759

The completed tree contains 476 tips and 391 alignment columns. Its artifacts,
input hash, tip identities, branch lengths and support ranges were checked.
All 952 focal-to-tip paths were independently recalculated by graph traversal;
maximum disagreement with the tree library was 1.33e-15.

For *Naganishia cerealis* (F610337), the original selected protein
KAJ9100180.1 overlaps an FCS EXCLUDE region. In this copy tree, an unrooted
five-tip split contains that sequence and *Candida albicans*, *C. dubliniensis*,
*C. maltosa* and *C. tropicalis*. The separating edge has reported SH-aLRT
support 100 and fitted length 0.2478307639.

Alternative protein KAJ9094313.1 has no recorded FCS coding overlap in the
existing candidate review, exact CDS translation, and 382 observed residues
across the 391 retained profile columns. Its unrooted three-tip split contains
*N. friedmannii* and *N. liquefaciens*, with reported SH-aLRT support 100 and
edge length 0.2857483549. A nested candidate–*N. friedmannii* split has support
67.4. The nearest tip by fitted path distance is *N. liquefaciens*; this is not
the same criterion as a sister-tip split.

This contrast supports prioritizing the alternative for orthology review. It
does not discriminate conclusively between mixed-source annotation, paralogy,
gene history and inference error. Support is SH-aLRT, not a bootstrap probability
or a probability of contamination. Absence of recorded FCS overlap is not a
cleanliness certificate; the *N. liquefaciens* assembly also has an unresolved
FCS publisher-checksum mismatch documented separately. No replacement, taxon
reassignment or new structure prediction has been made on this evidence.

Full paralog sampling, alignment/model sensitivity, rooting and reconciliation
remain necessary before accepting a replacement. The other copy trees remain
in progress. All focal split sides of size 2–8 and the ten nearest path-distance
neighbors are retained in the review, rather than selecting only the two
contrasting splits described here.

Reproduce the completed-tree review with:

```bash
python scripts/review_completed_fcs_copy_tree.py \
  --trees results/phylogeny/fcs-candidate-copy-trees-v1 \
  --inputs results/phylogeny/fcs-candidate-copy-tree-inputs-v1 \
  --marker 4764044at2759 \
  --output results/qc/fcs-copy-tree-review-4764044at2759-v1
```

Use a fresh output directory for re-execution. Full receipt and neighborhoods:
`metadata/fcs_copy_tree_review_4764044at2759.json`. Candidate sequence provenance:
`metadata/fcs_alternative_candidate_review.tsv`; retained alignment coverage:
`metadata/fcs_alternative_profile_alignment_review.tsv`.

## Second completed marker: 345792at2759

This tree contains 484 tips and 226 alignment columns. Its two alternative tips
and their selected counterparts yield 1,936 independently checked focal-to-tip
paths (maximum discrepancy 8.89e-16).

For *Acaulospora colombiana*, selected protein CAG8561170.1 overlaps FCS EXCLUDE.
It is on a four-tip unrooted split with three *Serendipita* entries (reported
SH-aLRT 83.7). Alternative CAG8780007.1 is instead nearest by tree path to
*Diversispora epigaea*, *Ambispora leptoticha* and *A. gerdemannii*. Its two-tip
split with *D. epigaea* has support 75.8. These neighborhoods differ, but support
for nearby alternative-copy splits is variable (including 35.2 for a five-tip
side); no fully resolved species placement is claimed.

For *Naganishia cerealis*, selected protein KAJ9101333.1 also overlaps EXCLUDE
and is nearest to *Candida* entries. Its two-tip split with *C. viswanathii* has
support zero and should not be treated as established. Alternative KAJ9091994.1
is nearest to *N. vaughanmartiniae*, *N. liquefaciens* and *N. friedmannii*. A
five-tip split contains the alternative, those three *Naganishia* entries and
*Filobasidium floriforme* (support 97.4). As in the first review, nearest path
distance is distinct from a sister-tip relationship.

These are two marker-level observations for N. cerealis and one for A. colombiana,
not independent proof of assembly contamination or candidate orthology. The
existing source-quality caveats, broader paralog/sequence checks, rooting and
reconciliation requirements remain. No candidate was promoted or predicted.
Reproduce using the command above with marker 345792at2759 and its own fresh
output directory. Receipt: metadata/fcs_copy_tree_review_345792at2759.json.

## Third completed marker: 4803900at2759

This tree contains 419 tips and 106 alignment columns. The selected and
alternative N. cerealis tips yield 838 independent path checks (maximum
discrepancy 1.78e-15). Selected protein KAJ9100835.1 overlaps FCS EXCLUDE
and lies on a four-tip split with C. albicans, C. dubliniensis and C. maltosa
(SH-aLRT 98.7). Alternative KAJ9112374.1 shares a two-tip split with
N. liquefaciens (SH-aLRT 100; separating branch length 1.0380).

The short alignment, long separating branch and unresolved N. liquefaciens
FCS checksum discrepancy require particular care. The observed contrast adds
a third marker for N. cerealis; these observations do not provide independent
proof of contamination or resolve candidate orthology. No replacement has
been accepted. All neighborhoods and provenance are retained in
metadata/fcs_copy_tree_review_4803900at2759.json. Reproduce with the same review
command using marker 4803900at2759 and a fresh output directory.
