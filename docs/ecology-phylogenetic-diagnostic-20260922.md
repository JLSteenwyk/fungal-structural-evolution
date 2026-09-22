# Ecological classification and minimum changes on species trees

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
