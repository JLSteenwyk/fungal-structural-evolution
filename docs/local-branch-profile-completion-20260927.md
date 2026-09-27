# Complete local codon branch profiles and unconstrained restart inventory

All 1,632 local MG94 cases completed eight fits each: one unconstrained
reoptimization and seven fixed target-branch parameter profiles at multipliers
0.1, 0.25, 0.5, 1, 2, 4 and 10. All 13,056 saved fits were reloaded in fresh
HyPhy processes for likelihood checks. The full audit checked 65,280 saved
artifact hashes and 224,752 fitted parameter values; maximum likelihood
readback discrepancy was 7.64e-11. Producer and auditor both terminated
successfully. Their source bindings and completed summary hashes were rechecked.

Ten cases have a constrained-grid likelihood exceeding the unconstrained result
by more than 1e-5. All ten are in Malassezia; the largest advantage is 6.9283
log-likelihood units for `Malassezia__5001413at2759__code1`. The next largest is
6.3895 for `Malassezia__4940884at2759__code1`. An unconstrained optimum of the same
model should be at least as good as an admissible constrained point. These
differences therefore identify optimization concerns; they are not likelihood
ratio evidence of selection or a biological explanation for the affected genus.

The complete next-stage inventory retains all eight audited parameter vectors
for every case, not just the ten flagged cases. It contains 13,056 distinct
within-case starts and 224,752 parameter values. Every vector was checked
against its saved profile model, original parameter names, branch identity and
unchanged non-parameter model text, then checked again after serialization.
The original free-model file, profile seed file, expected starting likelihood,
target parameter and all relevant hashes remain explicit. No seed was dropped
because its likelihood was lower or because its case lacked an initial flag.

Future unconstrained optimization must load the original free model and assign
the complete saved parameter vector. Loading a constrained profile export
without removing its target constraint would not be an unconstrained restart.
The starting likelihood must reproduce the audited profile value before
optimization, and every final model must be exported and checked in a fresh
process. No new optimization has been performed by the inventory stage.

Script: `scripts/prepare_local_mg94_multistarts.py`.
Starts: `results/cds/local-mg94-unconstrained-starts-20260927-v1`.
Full concern table: `grid_better_than_unconstrained.tsv` in that directory.
Completion/provenance: `metadata/local_branch_profiles_completed_20260927.json`
and `metadata/local_mg94_multistart_inventory_completed_20260927.json`.
Original full audit: `results/cds/local-mg94-branch-profile-audit-20260927-v1`.

The profiles concern the model's branch parameter t, while nuisance
exchangeabilities and omega change. Original synonymous-distance labels are
review metadata, not the profiled parameter. Finite grids and single starts do
not establish global optima, confidence intervals, absence of saturation or
eligibility for selection tests. Biological and contamination review remain
separate; zero recorded FCS action exposure is not proof of no contamination.
