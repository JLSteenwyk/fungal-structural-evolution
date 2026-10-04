# Direct predictor coordinates and phylogenetic controls

The completed control compares all **643** exact-complete-sequence
AlphaFold/ESMFold model pairs at all **12** original masks. These represent
673 marker/taxon links, 78 markers and 21 selected fungi. All **7,716** rows
remain explicit: **5,171** geometry comparisons and **2,545** insufficient-coverage
dispositions. No coordinate pair was rejected. The 71 markers in ready branch
controls are a different subset from the 78 with full-sequence model links.

## Methods and interpretation

Raw mmCIF checks validate complete polymer/all-atom/CA identity and confidence
against pinned models and native encodings. Both predictions must have valid
native features under the original six-residue confidence mask: pLDDT 0/70/90
and PAE unfiltered/5/10/15 Å. Comparison requires at least 50 residues and half
the full protein. Coverage exclusions remain missing, not zero displacement.

Measure proper-rotation CA superposition RMSD and distance MAE/RMS differences
for all residue pairs, spatial contacts, and sequence-local pairs. Spatial
contacts use the union of distances ≤15 Å in either model with original sequence
gap ≥3. Sequence-local pairs use original gaps 3–10; confidence gaps do not
renumber positions. Superposition excludes reflection; see
[Kabsch (1976)](https://doi.org/10.1107/S0567739476001873).

At pLDDT70/PAE10, **545** pairs have median CA RMSD **1.0378 Å**, all-pair distance
RMS difference **0.7064 Å**, spatial-contact difference **0.3838 Å**, and
sequence-local difference **0.2667 Å**. These are conditional predictor differences
for identical sequences. Whole-protein measures can reflect domain orientation;
neither experimental accuracy nor inherited evolutionary change follows.

The [PNG](figures/matched_predictor_coordinate_benchmark_20261004.png) and
[PDF](figures/matched_predictor_coordinate_benchmark_20261004.pdf) show state/local
geometry, local/whole-protein geometry, empirical distributions and all coverage
counts. Twelve threshold summaries and 96 descriptive rank correlations are
exported. Quantiles are not calibrated intervals; correlations have no p-values
or causal interpretation. Shared taxa/families/thresholds are not independent.

## Complete phylogenetic linkage

The join retains all **125 markers × 70 candidate tree views = 8,750 cases** from
the 526-taxon framework: 4,970 ready cases and 3,780 insufficient-observation
cases. All 4,523,750 original internal branch slots are accounted for:

| Original disposition | Branch slots |
|---|---:|
| Insufficient matched taxa | 1,954,260 |
| No observations on one side | 2,346,889 |
| Terminal projection | 50,544 |
| Unique internal projection | 13,235 |
| Shared internal projection | 158,822 |

The **31,290** observed internal paths include **13,235** paths mapping to one
original branch and **18,055** paths combining several original branches. Exact
original indices/units remain explicit; a merged path cannot locate change on
one original branch. Both sides retain taxon/model-pair identities and coordinate
coverage. Each path links its AA point estimate and six paired predictor
distributions (three models × two resampling modes). The 187,740 distribution
occurrences and 46,900 coordinate-link occurrences repeat existing controls;
they are not additional independent fits or replicates.

Whole-protein coordinate masks differ from marker-column branch-fit masks.
No cases were filtered by geometry. No geometric additivity, substitution-to-Å
conversion, time conversion or single-branch localization is accepted.
Candidate coalescent units remain distinct from substitution units and time.

## Verification, resources and reproduction

Independent raw-coordinate/scalar-mask/geometry replay checks every row using
`Rotation.align_vectors`, `pdist` and scalar sums without producer validation,
geometry or mask helpers. The standard Biopython mmCIF tokenizer is shared;
lexical independence is not claimed. Maximum geometry disagreement is
2.49e-14 Å, below predeclared 1e-10 arithmetic tolerances. Separate sorted
quantiles/tied ranks/scalar moments verify all summaries (maximum error
1.11e-16). A taxon-set reader independently reconstructs all original branch
slots and checks every case/path/source-data join. The figure was visually
inspected. These arithmetic tolerances do not change the older weighted guard.

All six original bounded waits, whole wrapper payloads, invocation-linked
manager start/end and **6,070** source/output bindings close. Each used two CPU
quotas, 16 GiB cgroup/12 GiB address-space caps, one BLAS thread and no swap.
Production/readback took about 169/175 seconds; aggregate six-stage child CPU
was about 373 seconds. No new fits, samplers, GPU predictions or charges occurred.
An initial software metric-label mismatch was corrected before full-stage
qualification/launch and retained in the qualification receipt. Identity,
proper rotation, reflection, scale, original-gap and invalid-input controls pass.

[Complete closure](../metadata/predictor_coordinate_benchmark_and_tree_links_completed_20261004_v1.json),
[coordinate reader](../metadata/overlap_coordinate_readback_20261004_v1.json),
[summary reader](../metadata/overlap_coordinate_summary_readback_20261004_v1.json), and
[path reader](../metadata/predictor_coordinate_phylogeny_link_readback_20261004_v1.json)
record exact checksums and execution identities. Raw datasets remain outside Git:

- `results/phylogeny/overlap-coordinate-benchmark-20261004-v1/all_model_pair_coordinates.tsv`.
- `results/phylogeny/overlap-coordinate-summary-20261004-v1/`: summaries/correlations.
- `results/phylogeny/predictor-coordinate-phylogeny-links-20261004-v1/`: all-case/all-path gzip JSONL.

Restore pinned raw models, encodings, original tree/marker inputs, point fits
and paired distributions at recorded paths; verify hashes. In a clean workspace
without immutable output roots, execute each command in its stage's
`metadata/*_execution_20261004_v1.json` using associated resources/wrapper.
Scripts refuse existing roots/receipts. Run software qualification, coordinate
production/readback, summary/readback, phylogenetic linkage/readback, then
`close_predictor_coordinate_benchmark_v1.py`. Close each original stage wait
with `record_software_original_wait_transport_v1.py` before its dependent reader.
Local custody does not establish public data availability; large-data release
remains a deliverable.

## Outstanding scientific requirements

The overlap has 21 fungi and no outgroups; it does not replace the
501-fungus/25-outgroup design. Full weighted fits remain zero and gated on
original full timing closure. Accepted species/gene/reconciliation/dating
frameworks, calibrated inference, broader predictor/experimental controls,
full atlas coverage and adequate ancestral posteriors remain required.
All eight aims remain incomplete. This benchmark does not repair historical
formatter/sampler errors or the SVD input-memory mutation.
