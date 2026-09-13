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
