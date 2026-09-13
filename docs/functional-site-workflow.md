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
