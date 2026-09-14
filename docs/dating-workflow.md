# Dating evidence and branch-duration requirements

The project currently reports sequence/structural divergence and conditional
model-relative site rates. Absolute change per unit time requires dated branch
durations with uncertainty, validated correspondence to each analysis tree and
a justified calibration model. Dating is an outstanding deliverable.

## Published chronogram inventory

The primary comparison is [Szánthó et al. (2025)](https://www.nature.com/articles/s41559-025-02851-z),
which explored fossil calibrations and relative time constraints across four
relaxed-clock scenarios. Its definition of Fungi differs from this project's
broad sampling convention; named clades must be matched by descendant membership,
not transferred by the label “Fungi”.

The existing [versioned Figshare archive](https://doi.org/10.6084/m9.figshare.28046594.v1)
was reverified against the publisher MD5. All archive files were inventoried with
uncompressed sizes and SHA-256 hashes. The four `*.timetree.tree.combined` files
each contain one annotated summary tree with 153 terminal and 305 total nodes.
Every annotated interval was matched to its node-age summary, and all 1,216
mean-age edge differences agreed with reported branch durations to within
4.55 × 10⁻¹³ in the published age units (Ma).

Each accompanying name map has 182 entries, of which 29 are absent from its
chronogram. These extra entries are explicitly marked and excluded from the
tip-overlap count. Each tree has 57 exact-name overlaps with the 526-entry project
manifest. Exact spelling is only a candidate correspondence: synonym, strain,
assembly and clade identity require review. Internal numeric labels are node
indices, not branch-support values.

These four files are mean chronograms with marginal uncertainty annotations.
They are not joint dated posterior samples. Marginal age intervals and pairwise
ordering probabilities cannot recover the joint distribution of parent–child
age differences. This inventory makes no claim that all possible posterior
artifacts from the authors have been located. No published dates or priors have
been assigned to project nodes.

```bash
python scripts/inventory_published_chronograms.py --output results/phylogeny/published-chronogram-inventory-v1
```

The receipt, all scenario counts and all name-map dispositions are versioned in
`metadata/published_chronogram_*`. The full archive-member inventory remains
outside Git with the existing 56.7 MB source archive.

## Calibration qualification and analysis requirements

Following [Parham et al. (2012)](https://doi.org/10.1093/sysbio/syr107), a proposed
fossil calibration needs a traceable specimen, evidence for phylogenetic
placement, reconciliation of relevant phylogenetic conflicts, explicit
stratigraphic context and a justified numerical age. A published estimate of a
node's age is not itself a primary fossil calibration.

Before a calibration is used here, record those fields together with:

- The exact project descendants defining the candidate node, whether the fossil
  constrains its stem or crown, and whether that node persists under the
  supported topology and contamination/taxon-identity sensitivities.
- The proposed minimum/maximum bounds or prior distribution, age uncertainty,
  original supporting sources and any alternative placements or superseded ages.
- Whether information is a primary fossil constraint, an inferred relative
  ordering, a secondary molecular age or an ecological scenario. Reused evidence
  must not be counted as independent calibration information.

The dating analysis should compare defensible calibration sets and clock/model
choices and examine the prior-only distribution as well as the posterior.
Adequate sampling and convergence must be checked for the actual analysis.
Structural divergence must not be used both to set branch duration and to test
for structural acceleration on that branch.

Branch-level time-normalized comparisons must propagate matched parent/child
ages from joint samples or a justified alternative that retains their covariance.
Dividing a point estimate by a mean branch duration does not propagate dating
uncertainty. Mapping durations onto reconciled gene trees also requires explicit
handling of duplication/loss and different descendant sets. Until these
requirements are met, figures and tables must continue to label quantities as
relative divergence or conditional rate multipliers, without units of change
per million years.

## Published calibration candidate catalogue

Downloaded and hashed the publisher's
[Supplementary Information sections 1–6](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-025-02851-z/MediaObjects/41559_2025_2851_MOESM4_ESM.pdf).
Section 1 contains 27 calibration headings in numbered groups 1–24, with A/B
entries for groups 7, 22 and 24. The catalogue records 21 published minimum
bounds and six soft maximum bounds, their published ages, stem/crown/root labels,
anchor names, PDF page and exact-name project matches. A/B entries are kept
separate; they must not automatically be treated as constraints on identical
nodes. Published geological ages have not been updated or approved for use here.

A second PDF extraction engine (PyMuPDF) independently located each heading on
its reported page and checked all 27 numeric ages and bound types against the
pdftotext catalogue. Warnings from the full pdftotext extraction are preserved.
This verifies factual transcription, not specimen identity, geological age
justification or the calibration's phylogenetic placement.

Four entries have both anchor names exactly represented in the current project
manifest. Lack of an exact pair does not itself disqualify a calibration:
representative descendants may define the same clade after explicit node mapping.
Conversely, matching both names does not resolve the meaning of a published stem
versus crown label, outgroup rooting, accession identity or topology uncertainty.
No calibration has yet been accepted and no clock inference has been launched.

```bash
python scripts/catalog_published_fungal_calibrations.py --output results/phylogeny/published-calibration-candidates-v1
```

Reproduction requires the PDF at the path in
`metadata/published_calibration_candidates_receipt.json`, downloaded from its
recorded publisher URL and checked against the recorded SHA-256. Structured
candidate facts and independent extraction evidence are versioned in
`metadata/published_calibration_candidates.tsv` and
`metadata/published_calibration_heading_readback.json`. Full source PDFs/text
remain outside Git. Next steps are specimen/source review and descendant-based
mapping of each candidate onto qualified project topologies.
