# Candidate functional-site correspondence

Pfam 38.2 active-site annotations have been projected to the full marker-protein
set. The active-site file uses HMM positions, not full UniProt sequence
coordinates; Pfam describes this mapping scheme in its [release explanation](https://xfam.wordpress.com/2017/03/08/pfam-31-0-is-released/).
The pinned file, release metadata, source profile HMMs and complete significant
hit annotations are checked before projection. The parser retains source
accessions, full reference patterns and selenocysteine annotations.

All 10,570 significant hits belonging to 173 annotated families were processed,
covering 8,961 distinct marker protein sequences. Each hit's aligned protein
segment was realigned to its own versioned Pfam HMM using HMMER 3.4. Stockholm
reference annotations identify profile match columns; insertion and deletion
positions remain explicit. The entire ungapped segment must match its source
protein exactly, and each projected residue is checked against its full protein
coordinate. Tests cover insertions, deletions, offsets, changed sequences,
out-of-profile positions and the pinned annotation format.

A conservative candidate flag requires every annotated residue in one reference
pattern to be conserved, at least 50% coverage of the original profile hit and
no alignment overlap with another retained gathering-threshold hit. The 50%
rule is an operational screen, not calibrated functional confidence. Family and
Domain annotation types remain distinct in the full table. All overlapping,
substituted and unaligned rows are retained. Unconserved patterns do not prove
loss of function, and conserved patterns do not prove catalytic activity.

The output contains 19,429 reference-pattern/site rows, including 15,704 mapped
residue coordinates. Independent readback checked all those coordinates and
13,655 whole-pattern statuses. There are 5,258 candidate rows, reducing to
**5,251 distinct protein positions in 3,503 proteins across 94 Pfam families**.
Repeated reference annotations are not independent site observations.

These are candidate functional correspondences. They are not experimental
validation, resolved orthology, structural-confidence qualification or tests
of selection. Joining the positions to prediction-specific structures and
sequence/structural phylogenies is the next step; uncertainty and multiple
testing remain required before functional evolutionary claims.

```bash
python scripts/map_pfam_active_sites.py \
  --output results/functional_sites/pfam-marker-v1
python -m unittest discover -s tests -p test_pfam_active_sites.py
```

The complete row-level table, profile alignments and provenance are in the
immutable result directory. The deduplicated candidate index, receipt and
independent readback are versioned as `metadata/pfam_active_site_*`. This first
projection covers the complete marker annotation set; whole-proteome extension
awaits the complete domain searches.

## Structural and paired-alignment joins completed

All projected correspondences, not only conserved candidates, were joined
separately to the expanded AlphaFold and frozen local ESMFold snapshots.
Reference-pattern duplicates were aggregated to 16,849 hit/profile-site
identities while retaining expected residues, reference accessions and all
pattern statuses. Protein links yield the same 17,105 taxon–marker–site rows
for each source, including gapped, nonconserved and ambiguous correspondences.
Keeping these rows avoids conditioning future analyses of functional-site
change entirely on present-day residue conservation.

The joins distinguish model availability, a mapped protein coordinate, valid
native structural features, joint six-residue confidence and actual observation
in the emitted paired phylogenetic alignment. Invalid native states remain
unobserved. A site can have structural confidence but lie outside the retained
marker matrix or an eligible paired alignment; these cases remain explicit.

AlphaFold contributes 2,713 observed rows across 244 taxa and 21 markers,
including 1,119 conserved candidates. ESMFold contributes 1,189 observed rows
across 191 taxa and 13 markers, including 139 conserved candidates. The union
contains 3,902 observations across 399 taxa, with no shared observed
site/taxon/marker identities between these acquisition snapshots. This is
coverage of homologous correspondences, not 3,902 independent evolutionary events.

Independent readback verified the complete identical row universe, unique
identities, joint confidence and all observed amino-acid/3Di characters at their
reported paired-alignment columns. Three tests cover retention of nonconserved
and gapped sites, reference-pattern deduplication and rejection of conflicting
projections. Full result tables retain unavailable and excluded sites; versioned
`metadata/{gdm,esmfold}_functional_sites_paired.tsv` files contain the observed
subsets, and receipts pin the full outputs. Sources remain separate for inference.

```bash
python scripts/link_functional_sites_structures.py \
  --functional-sites results/functional_sites/pfam-marker-v1 \
  --snapshot results/structural_markers/gdm-expanded-v1 \
  --encodings results/structural_alphabet/audited-gdm-expanded-v1 \
  --paired results/phylogeny/paired-inputs-gdm-expanded-v1 \
  --source-label AlphaFold --output results/functional_sites/gdm-linked-v1
python scripts/link_functional_sites_structures.py \
  --functional-sites results/functional_sites/pfam-marker-v1 \
  --snapshot results/structural_markers/esmfold-partial-v1 \
  --encodings results/structural_alphabet/audited-esmfold-partial-v1 \
  --paired results/phylogeny/paired-inputs-esmfold-partial-v1 \
  --source-label ESMFold --output results/functional_sites/esmfold-linked-v1
python -m unittest discover -s tests -p test_functional_site_join.py
```

Branch localization, ancestral states, tests against matched background sites,
functional enrichment and experimental interpretation remain pending. The
joins provide checked coordinates and confidence; they do not themselves
establish structural acceleration, positive selection or catalytic activity.


## Functional annotations in the local site-evolution frame

Completed `results/phylogeny/site-evolution-functions-esmfold-v1` retains all
16,571 original sites across 72 local-source markers and adds the existing
observed functional correspondences. The join verifies common paired-input
provenance, model/protein identity, residue and matrix/paired column positions,
and observed AA/3Di states against the accessibility projection.

The 1,189 annotation rows span 1,066 distinct taxon-site observations at 28
paired sites across 13 markers. Multiple profile correspondences can refer to
one observed residue; per-site taxon counts therefore use sets. Conserved
candidates contribute 139 annotation rows at 10 sites across seven markers.
Candidate and other-correspondence taxon sets may overlap. The unannotated
sites remain explicitly `no_observed_correspondence`; this is not evidence of
nonfunction. Fractions use all observed taxa, not an assumed annotation
sensitivity or a count of independent evolutionary events.

All inherited frame and annotation fields passed full readback. Every appended
coordinate ASA/confidence/context field, including repeated taxon-site rows,
matched its source. All taxon sets, fractions, Pfam accession unions, status
fields and marker summaries were reconstructed. These annotations enable
later focused evolutionary case studies; no functional enrichment, branch
effect, selection test or experimental functional validation has been performed.
The limited current coverage does not support generalizing to fungal functional
sites as a whole.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/annotate_site_evolution_functions.py --frame results/phylogeny/site-evolution-frame-esmfold-v1 --functional results/functional_sites/esmfold-linked-v1 --projection results/structural_annotations/paired-accessibility-esmfold-v1 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --output results/phylogeny/site-evolution-functions-esmfold-v1
```

Full tables remain outside Git; receipt, readback and marker coverage are in
`metadata/esmfold_site_evolution_functions_*`.

## Complete ESMFold-cohort functional join running — 23 September 2026

The earlier local-source functional join remains an acquisition checkpoint.
A new full-cohort join now uses all 25,322 completed ESMFold models and the
122-marker paired input set. It preserves the complete original functional
correspondence universe, including nonconserved and gapped sites. It does not
pool prediction sources or reinterpret conserved patterns as proven activity.

`scripts/link_completed_esmfold_functional_sites.py` preserves the existing
join's coordinate, sequence, native-validity, confidence and paired-character
checks while requiring a completed qualified-encoding union and the full
paired-array readback. Every requested mapping/encoding/paired receipt must
be covered by that readback and retain its exact hash. The joint pLDDT/PAE
qualification stage is required. The original producer remains unchanged for
reproducibility of the older checkpoint.

The three existing tests for retained nonconserved/gapped rows, reference
pattern deduplication and conflicting projections passed against the new
script's aggregation function. Full production readback remains pending.
The live unit is `fungal-completed-esmfold-functional-sites-20260923.service`;
its launch record and prelaunch source/resource plan are in
`metadata/completed_esmfold_functional_site_join_{launch,plan}.json`.
One CPU and 8 GiB RAM are allowed, without swap, GPUs or new charges. Output
planning is 1 GiB and runtime planning is an uncalibrated 0.02–1 hour.

```bash
python scripts/link_completed_esmfold_functional_sites.py \
  --functional-sites results/functional_sites/pfam-marker-v1 \
  --snapshot results/structural_markers/esmfold-all-completed-20260922-v1 \
  --encodings results/structural_alphabet/audited-esmfold-all-completed-20260922-v1 \
  --paired results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1 \
  --paired-readback metadata/esmfold_all_completed_paired_inputs_completed_readback.json \
  --source-label ESMFold \
  --output results/functional_sites/esmfold-all-completed-linked-20260923-v1
```

This refresh is running, not yet a completed functional result. Branch/site
localization, matched background tests, ancestral uncertainty and experimental
interpretation remain outstanding.

## Complete ESMFold functional join and full readback passed

The full-cohort job finished with exit status zero. Independent reconstruction
now passed **every field of all 17,105 taxon–marker–site rows**, retaining
16,849 unique profile-site correspondences. The completed ESMFold cohort
provides **5,576 paired-observed rows across 283 taxa and 21 markers**:
2,240 conserved-candidate rows and 3,336 other correspondence rows. These
are annotation rows, not necessarily distinct protein residues or independent
evolutionary events. The earlier 1,189-row local-source coverage remains an
older acquisition checkpoint and should not be presented as current coverage.

The complete output retains 7,783 rows without a model, 282 modeled sites
outside the marker matrix, 5,768 mapped rows and 3,272 gaps in the functional
profile alignment. Of the mapped rows, 5,576 are observed in the qualified
paired alignments. Those availability categories are not functional-loss
classes and do not indicate whether a missing model exists in another source.

`scripts/readback_completed_functional_sites.py` reconstructs the reference
pattern aggregation and complete row universe from original projections and
protein links, then checks source model identities, residue/matrix positions,
all emitted confidence/state fields and both presence and absence in paired
FASTA arrays. It imports no producer join functions, while sharing existing
source projections, qualified arrays, parsers and the checksum helper. It does
not repeat HMM alignment, coordinate reconstruction or PAE qualification.

The completion and readback receipts are archived as
`metadata/completed_esmfold_functional_site_join_{receipt,readback}.json`.
The complete result table remains at
`results/functional_sites/esmfold-all-completed-linked-20260923-v1/site_structure_links.tsv`.
Reproduce its readback with:

```bash
python scripts/readback_completed_functional_sites.py \
  --result results/functional_sites/esmfold-all-completed-linked-20260923-v1
```

This supersedes the running/readback-pending status immediately above.
Functional-site evolutionary effects, matched background tests, ancestral
uncertainty and experimental interpretation remain unfinished.

## Recovered AlphaFold functional correspondences — September 27

The updated join retains the same complete universe of 17,105 annotation rows
and 16,849 unique profile-site identities. It uses the recovered AlphaFold
mapping/encoding union and current 125-marker paired inputs, with the completed
paired-array readback required. All projected protein/model identities,
coordinates, native states, confidence fields, and presence or absence in the
paired matrices passed the independent row-by-row checker.

There are **6,444 observed correspondence rows across 314 taxa and 23 markers**,
including 2,555 conserved-candidate and 3,889 other correspondence rows. This
supersedes the earlier 2,713-row AlphaFold snapshot for coverage reporting.
The full table retains 6,723 rows without a model, 356 modeled rows outside the
marker matrix, 6,754 mapped rows, and 3,272 gaps in the functional profile
alignment. None of these categories establishes functional loss.

The producer is `scripts/link_recovered_afdb_functional_sites.py`, a separate
version of the completed-cohort join requiring the recovered AlphaFold encoding
receipt. The original ESMFold script and results remain unchanged. Arguments,
source pins and resource estimates are in
`metadata/recovered_afdb_functional_site_join_plan_20260927.json`; exact process
identities are in the adjacent join and readback launch records. Output is
`results/functional_sites/afdb-recovered-linked-20260927-v1`. Recheck with:

```bash
python scripts/readback_completed_functional_sites.py \
  --result results/functional_sites/afdb-recovered-linked-20260927-v1
```

## Functional correspondence and accessibility integration

`scripts/annotate_recovered_site_functions.py` joins the observed annotation
rows to the audited raw accessibility projection and site parsimony/exposure
summaries. Both full source readbacks are mandatory. The output retains all
47,529 sites across 125 markers and all original fields, adding annotation
counts, distinct taxon counts, explicit missing-annotation status and Pfam IDs.
The annotation ledger adds each observed residue's isolated-chain ASA and
confidence. Candidate and other-correspondence taxon sets can overlap; they
are not added together as distinct observations.

Output is
`results/functional_sites/site-parsimony-exposure-functions-afdb-recovered-20260927-v1`.
The complete CLI and pinned inputs are in
`metadata/recovered_afdb_functional_exposure_plan_20260927.json`. The production
join has completed: correspondences occur at 50 paired sites, including 18
sites with conserved candidates. Full independent verification passed
using `scripts/readback_recovered_functional_exposure.py --plan
metadata/recovered_afdb_functional_exposure_plan_20260927.json`.

The readback independently compares every inherited site value, reconstructs
all observed annotation rows and coordinate joins, and checks every added
per-site count/fraction and marker summary. These stages each use one CPU,
8–16 GiB memory, no swap and at most 1 GiB planned output; 0.05–2 hours are
reserved per join/readback, with no GPU work or paid infrastructure.

These are sparse projected functional correspondences, not validated catalytic
sites or binding pockets. The 6,444 rows are not independent evolutionary
events. Sites without observed annotations remain unknown, not nonfunctional.
Parsimony counts are minimum changes, not rates, branch assignments or selection.
Matched-background tests, rate integration, ancestral uncertainty, functional
enrichment and experimental interpretation remain open.

The complete integration readback passed for all 47,529 sites and 6,444 ledger
rows (6,190 distinct annotated taxon–site observations). All four producer/audit
services are terminal with exit status zero. Source pins and declared artifacts
were rechecked after completion. Archived receipts and terminal-state evidence
are in `metadata/recovered_afdb_functional_exposure_completed_20260927.json`;
the source results remain outside Git. This supersedes any running status for
these two joins, not the outstanding functional evolutionary analyses.

## Complete predictor comparison at annotated residues

The recovered AlphaFold and completed ESMFold functional joins were compared
across the identical full 17,105-row annotation universe. Shared protein,
reference-pattern and residue fields must agree. Predictor-specific availability,
models, confidence, matrix/paired columns and native states remain separate.
A complete dataframe join independently checked every exported source field;
all coverage/state classifications and exact-coordinate aggregation were also
read back. Reproduce with:

```bash
python scripts/compare_functional_prediction_sources.py \
  --output results/functional_sites/prediction-source-comparison-NEW
```

The completed output is `results/functional_sites/prediction-source-comparison-20260927-v1`.
The full annotation table retains gapped positions. The separate residue table
collapses repeated profile correspondences at the same marker/taxon/protein/
coordinate, retaining whether any annotation is a conserved candidate and
whether any is another correspondence; these flags can overlap.

| Paired observation coverage | Annotation rows | Distinct protein coordinates |
|---|---:|---:|
| Both predictors | 158 | 150 |
| AlphaFold only | 6,286 | 6,040 |
| ESMFold only | 5,418 | 5,199 |
| Neither | 5,243 | 1,874 |

The coordinate table contains 13,263 rows; annotations without a protein
coordinate remain in the full table rather than being assigned artificial
residue identities. The 150 jointly observed coordinates span 14 taxa,
14 markers and 95 taxon/protein identities. Their structural-alphabet states
agree at 127 coordinates and differ at 23. Of 59 jointly observed coordinates
with any conserved-candidate annotation, 45 agree and 14 differ. These are
conditional counts, not independent observations or an enrichment test.

Fourteen of the 23 differing coordinates belong to marker `5005892at2759`,
whose differing annotation rows map to Pfam `PF01163.29` (RIO1). This concentration
is a predictor-sensitivity follow-up, not evidence of a lineage-specific change.
The exhaustive 23-coordinate review list is archived as
`metadata/functional_prediction_discordant_residues_20260927.tsv`; denominator
and source proofs are in
`metadata/functional_prediction_source_comparison_completed_20260927.json`.

Missing paired coverage is never treated as a state difference. Structural-
alphabet states depend on nonlocal residue context; agreement does not establish
accuracy and disagreement does not identify which predictor is correct.
Prediction circularity, correlated homologous sites and sparse overlap prevent
these counts from providing broad independent functional validation. Direct
coordinate/context comparisons and experimental evidence remain needed.

## Full comparison of native residue contexts

All 150 jointly observed functional coordinates were examined in both
predictors' qualified arrays, including the agreeing states. The comparison
loaded 186 distinct qualified model/source entries. Exact full-sequence identity,
focal residue/state, native partner index and the six role positions
`i−1, i, i+1, j−1, j, j+1` were checked. Every position stays in bounds;
minimum context pLDDT was recomputed from the C-alpha confidence array and
checked against both the qualified encoding and the functional join. Every
context retains pLDDT ≥70 and audited context PAE ≤10. PAE maxima were read from
previously qualified arrays, not independently recomputed in this comparison.

| State comparison | Same partner position | Different partner position |
|---|---:|---:|
| Same state | 118 | 9 |
| Different state | 11 | 12 |

Thus 12 of the 23 state differences accompany a partner change, while 11 do not.
The nine same-state/different-partner observations remain in the analysis.
A changed partner is not proof that partner selection caused a state change;
the same partner does not mean that coordinates or descriptors agree. These
are conditional observations in a small, correlated subset, not error rates.
Direct coordinate/descriptor comparisons remain necessary before interpreting
the affected residues as evolutionary or functional changes.

Reproduce with `python scripts/compare_functional_structural_contexts.py --output
results/functional_sites/prediction-context-comparison-NEW`. Completed output is
`results/functional_sites/prediction-context-comparison-20260927-v1`. All exported
contexts were reloaded and compared with their qualified arrays; this is a
serialization/source check, not fresh native encoding or coordinate validation.
The complete 150-row table and provenance receipt are versioned as
`metadata/functional_prediction_contexts_20260927.tsv` and
`metadata/functional_prediction_context_receipt_20260927.json`. Both source
models, six-position contexts, contextual amino acids and confidence summaries
are retained for every row. No new predictions or GPU work were performed.
