# Ecological classification and minimum changes on species trees

## Sample-linked extension to 32 species

The new `metadata/species_ecology_evidence_sample_linked_20260922.tsv`
preserves the original 26 statements and adds Amanita rubescens, Gautieria
morchelliformis, Hysterangium stoloniferum, Ramaria rubella and both sampled
Thelephora species. Each added classification is directly reported in
Miyauchi Dataset 2 and linked to the selected assembly by Dataset 1's
BioSample and WGS project. `scripts/curate_sample_linked_ecology.py` records
the curated choices and checks these source/identity joins. Five additions
are ectomycorrhizal; Ramaria is source-classified saprotrophic. Its disagreement
with the genus-level candidate and the source's R. acris naming remain
explicit. This is a source-specific coding, not a resolved taxonomy claim.

Rerunning the unchanged diagnostic with
`metadata/ecology_sample_linked_parsimony_plan.json` gives **seven minimum
changes** on both ML trees and all 2,000 bootstrap trees. The binary coding
has 24 ECM and seven asymbiotic/saprotrophic classifications; 495 tips remain
unknown, including the orchid symbiont. There are 2,064 output rows. The
independent graph-cut readback passes all 64 ML/leave-one-out checks.
Treating Ramaria's classification as unknown reduces both ML scores to six.
Thus the new score depends on its uncertain source-specific assignment;
it must not be presented as a robust count of seven independent origins.

The new output is
`results/ecology/sample-linked-classification-parsimony-20260922-v1`, with
versioned `ecology_sample_linked_parsimony_{receipt,readback}.json` records.
The frozen 26-species analysis below remains available as a sensitivity
comparison. Neither analysis establishes transition locations, ancestral
direction, ecological effects or sufficient independent replication.

## Earlier 26-species diagnostic

The 26-species evidence table now has an explicit phylogenetic diagnostic on
both completed profile-alignment C20-PMSF trees. The binary coding retains
19 published ectomycorrhizal classifications and six source-classified
asymbiotic or saprotrophic taxa. All other tips are unknown, including the
orchid symbiont Tulasnella: orchid association alone does not establish the
absence of ectomycorrhiza. In total, 501 of 526 tips remain unconstrained.
The source classifications do not prove that ecological capabilities are
exclusive or that the selected isolates were experimentally validated.

A symmetric, unit-cost Sankoff calculation requires a minimum of **six
total changes** on each maximum-likelihood tree. Every saved ultrafast
bootstrap tree also has a score of six (1,000 per guide, 2,000 total).
This is a conditional minimum-change diagnostic, not six independently
verified origins. Unknown tips can take either state; neither the root
state nor the gain/loss direction is specified. The bootstrap result
describes sampled topology variation on one sequence alignment, not
uncertainty in the trait classifications, transition model or alignment.

Leaving any one of the following classifications unknown lowers the
minimum to five on the profile-guide tree: Amanita inopinata, A. thiersii,
Botryobasidium botryosum, Tuber melanosporum, Cenococcum geophilum, Glonium
stellatum, Lepidopterella palustris, or Sphaerobolus stellatus. The full
52-row ML/leave-one-out sensitivity grid is independently verified by
graph minimum cuts. These influential classifications prioritize evidence
review; they do not justify excluding taxa to select a preferred result.

The dynamic program passed exhaustive binary assignments for 243 small-tree
cases, including unknown tips and a four-way polytomy. The independent
graph-cut calculation verifies both full ML scores and every leave-one-out
score. Bootstrap scores were computed and their tip identities and counts
checked by the primary program; they were not independently recomputed.

Reproduce the primary calculation into a fresh output directory specified
by a new plan if the original output already exists:

```bash
python scripts/map_ecology_parsimony_sensitivity.py \
  --plan metadata/ecology_classification_parsimony_plan.json
python scripts/readback_ecology_parsimony_min_cut.py \
  --plan metadata/ecology_classification_parsimony_plan.json \
  --output metadata/ecology_classification_parsimony_readback.json
```

The complete 2,052-row table is outside Git at
`results/ecology/classification-parsimony-20260922-v1/minimum_changes.tsv`;
its checksum and summaries are versioned in
`metadata/ecology_classification_parsimony_receipt.json`. Plans pin source
trees, bootstrap files, their prior readbacks, the sampling manifest and
the expanded ecological evidence. No branch is yet accepted as an
independent ecological transition. Further curation, ancestral-state model
sensitivity, completed crossed trees and matched structural coverage are
required before association testing.

The expanded confidence-qualified ESMFold inputs now have a separate
[full ecological overlap diagnostic](qualified-ecology-coverage-20260922.md).
It intersects actual alignment observations across the 32 curated species;
28 have eligible ESMFold markers. This coverage check does not change the
parsimony results or establish independent ecological transitions.
