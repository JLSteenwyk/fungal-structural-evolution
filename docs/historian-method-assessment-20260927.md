# Historian assessment for joint ancestral indel histories

The independent SIC estimates failed compatibility constraints in 614 of 918
candidate node/model combinations. Conditioning independent gap weights on
compatibility caused large probability changes. A coherent evolutionary history
model therefore needs assessment before ancestral sequences are qualified.

Historian is a candidate, not an accepted replacement. The pinned author source
is commit `359be7ef8e71a04c90aaffe44b0ddf6901a12ff0` from
[evoldoers/historian](https://github.com/evoldoers/historian).
Its README identifies the methods papers by
[Holmes (2017)](https://doi.org/10.1093/bioinformatics/btw791) and
[Westesson et al. (2012)](https://doi.org/10.1371/journal.pone.0034572).
This assessment is based on the pinned implementation and its documentation;
it does not claim an independent replication of their published benchmarks.

## Build and upstream tests

Source lives under `data/software_audits/historian-20260927/source`.
Clone the author repository there and check out the commit above before running
the build script. Dependencies are versioned Ubuntu GSL packages downloaded
with `apt-get download` and unpacked locally with `dpkg-deb -x`; no system
packages are installed. Exact package versions, compiler, commands, binary and
package hashes are retained in the build plan and receipt.

The original v1 build failed at linking because upstream places shared
libraries before objects. `build_historian_local_v2.py` adds
`-Wl,--no-as-needed` to retain these libraries, leaving source unchanged.
The successful build is `results/software/historian-build-20260927-v2`.
`test_historian_local.py` runs upstream `make -k test` with the same flags;
the suite exited zero. Logs and hashes are retained under
`results/software/historian-tests-20260927-v1`.
This establishes local build/test functionality, not project-scale capacity.

## Implementation findings that affect interpretation

- Branch transducers approximate indel evolution by excluding overlapping
  indel events within a branch (README, method section).
- The algorithm retains reduced ancestral profiles, so a profile state cap,
  sampling count or threshold, seed and alignment band must be explicit.
  Memory-derived defaults must not silently determine scientific settings.
- `recon.cpp` expands an infeasible guide band, ultimately removing the guide
  constraint. Even `-band 0` is therefore not a guarantee of fixed alignment.
  Logs must be checked for this behavior in every run.
- Reconstruction selects a best root alignment path. `-ancprob` then computes
  residue probabilities conditional on that reconstruction; it does not supply
  a full joint posterior over indel histories. The exported residue probabilities
  omit states below 0.01 (`sumprod.h`); absent entries are not zero probability.
- The MCMC implementation requires ultrametric trees. It is not directly
  applicable to the current gene trees. Artificially making them ultrametric
  would change the inference problem.
- `tree.cpp` requires bifurcating trees and parses branch lengths using a
  minimum of 1e-9 (`tree.h`). Root placement remains an assumption for joint
  insertion/deletion histories; reversible substitution fitting did not identify
  the suppressed degree-two root or its edge split.

## Complete project input audit

`audit_historian_project_inputs.py` checks all 78 whole/domain alignments in
13 families, preserving up to 622 proteins per alignment. It verifies exact
tip identities, alignment widths and pinned tree provenance. Results are in
`metadata/historian_project_input_audit_20260927.json`.

Six alignments of OG0000230 contain the same four-child node, n1790. Historian
would reject these inputs. Twenty-four alignments have nonroot branches below
the parser minimum; the audit records every affected node and original value.
No inputs were altered or removed, and no project inference was launched.

Before full reconstruction, evaluate a justified treatment of the polytomy
(including sensitivity across resolutions if appropriate), effects of the
branch floor, root placement, full-family memory/runtime and guide/profile
approximations. Model fitting and complete uncertainty propagation remain
required. Passing software tests alone does not establish ancestral reliability.

## Full sensitivity grid launched

Further inspection shows that the four children of n1790 are terminal proteins
with identical aligned sequences in all six affected inputs, and all four edges
have length zero. We retain every copy and enumerate all 15 rooted binary
resolutions. A separate enumeration by the balanced and caterpillar four-leaf
shapes verifies that the resolution set is exhaustive and duplicate-free.
The inserted internal nodes are explicitly artificial and are not biological
ancestor candidates.

`prepare_historian_sensitivity_inputs.py` generates derived trees at minimum
edge lengths 1e-9 and 1e-7. It retains source roots and checks all original nodes,
all tips, contraction of artificial nodes, every edge length and every change
to root-to-tip distance after serialization. The 108 tree files support 324
runs covering all 78 original alignments. Original files remain unchanged.
Derived trees and their complete checksummed mapping are stored under
`results/ancestral/historian-sensitivity-inputs-20260927-v1`.

`run_historian_capacity_grid.py` executes the full grid using fixed diagnostic
parameters: LG, four gamma categories with shape 1, insertion/deletion rate
0.01, expected gap length 3, guide band 20, 100 profile samples, maximum
20,000 profile states, seed 20260927 and no iterative refinement. These are
computational checks, not family-fitted evolutionary estimates. JSON outputs
contain reconstructed sequences, not a claimed posterior ensemble. Unknown
extant X positions are explicitly recorded if imputed; every known residue
and every tip length must be preserved. Output topology, root and branch
serialization are checked. Guide-band relaxation messages are retained.

Plan: `metadata/historian_capacity_plan_20260927.json`. Limits: two CPU workers,
32 GiB aggregate memory and no swap, with 12 GiB sampled RSS and 30 minutes
per process. A process that exceeds a per-job limit is terminated as a process
group and retained as an unresolved capacity result, not dropped or scored as
successful. Planning allowance is 1–82 hours and 10 GiB, with no paid resources
or GPUs. Per-job receipts make successful and unsuccessful dispositions
reviewable; overall numerical and biological validation remain downstream.

Root placement, model fitting, approximation sensitivity and joint history
uncertainty remain necessary even if every capacity check passes. Comparing
these equally enumerated resolutions is a sensitivity assessment, not a
posterior weighting of tree histories.

## Additional uncertainty limitation from source inspection

`ForwardMatrix::bestTrace` explicitly describes its traceback as not quite
Viterbi (`src/forward.h`). It selects predecessors using forward sums. In
`makeProfile`, an effective transition sums probabilities over eliminated
paths but retains one representative `bestAlignPath` for sequence output
(`src/forward.cpp`). Consequently, changing only the root traceback to random
sampling would not recover all eliminated descendant alignment histories.
Calling such a modification a joint-history posterior sampler would be
unjustified. No source modification or new sampler was made.

A second independent output snapshot checks 70 finished runs, leaving 254
pending. All available checks pass; all 1,532 one-factor candidate comparisons
have zero ungapped edit distance. This includes available star resolutions,
but is neither completion of the grid nor validation of posterior uncertainty.

## Observed memory failures and unchanged-model retries

The two OG0000972 whole-MAFFT runs (622 proteins; both edge floors) exceeded
the 12 GiB sampled RSS cap after approximately 12.6 and 11.9 minutes, with
recorded peaks of 12.02 and 12.29 GiB. They were killed by the declared resource
policy and remain failed capacity dispositions in the original grid.

`run_historian_memory_retries.py` retries these two frozen jobs sequentially at
64 GiB per-process sampled RSS and a two-hour timeout, inside an 80 GiB/no-swap
service with one CPU. Planning allowance: 0.5–4 hours and5 GiB, no paid resources.
Every sequence, tree, seed, rate, band, sample count and profile-state limit is
unchanged. The retry manifest pins the original failure receipts. Other initial
grid jobs continue; later failures require their own explicit disposition.
`audit_historian_memory_retries.py` provides the same independent output checks
for this separate retry collection. No retry completion has yet been claimed.

The paired whole-FAMSA runs subsequently also exceeded the 12 GiB cap, after
12.9 and 13.2 minutes (recorded peaks 12.05 and 12.03 GiB). Their original
failed receipts and input hashes are pinned in a separate FAMSA retry plan.
`run_historian_famsa_memory_retries.py` now runs both jobs sequentially with
the same 64 GiB RSS/two-hour limits and 80 GiB/no-swap service ceiling.
Planning allowance is 0.5–4 hours, 5 GiB output, and no paid resources. The
622 proteins, trees and scientific options match the original jobs exactly.
The existing MAFFT retry batch and original capacity grid are unchanged.
`audit_historian_famsa_memory_retries.py` is prepared for independent readback;
completion and successful reconstruction remain unproven until outputs pass.

The first whole-MAFFT retry (floor 1e-9) finished in 1,000 seconds with a
16.82 GiB peak sampled RSS. Independent readback passed for all 622 tips,
the output tree and all four candidate descendant sets. The paired MAFFT
retry remains pending, so no paired sensitivity comparison is yet available.
This is computational feasibility and output integrity, not a fitted or
converged ancestral posterior. Both frozen retry collections now have queued
final audits that require the exact producer identity, successful terminal
service state, all expected job receipts and artifact hashes. Unsuccessful
job outcomes remain explicit even if the batch driver exits successfully.
The first partial readback is recorded in
`metadata/historian_memory_retry_partial_readback_20260927.json`.

## Completed MAFFT retries reveal branch-floor sensitivity

Both 622-protein MAFFT retries and their independent readback are terminal
with successful exit status. Runs took 1,000 and 920 seconds and peaked at
16.82 and 16.61 GiB sampled RSS. All tree, extant-sequence and candidate
mapping checks pass. Original 12 GiB capacity failures remain recorded.

Changing the minimum branch length from 1e-9 to 1e-7 changes every candidate
output. Non-root levels 0, 1 and 2 have ungapped edit distances 10, 15 and 15;
lengths change 467→467, 461→467 and 460→467 respectively. The assumed-root
candidate changes 460→467 with edit distance 15. These are computational
sensitivity measurements, not inferred biological event counts. Same-length
edit distance does not by itself establish substitution counts.

This differs from the earlier small-family snapshots with no observed edits.
The approximate profile/traceback reconstruction is not robust to this
numerical floor choice for the largest family. Both alternatives must remain
visible; neither is selected as a qualified ancestor on these results. The
FAMSA retry collection and joint-history sampling diagnostics remain separate.
Evidence: `metadata/historian_memory_retry_audit_completed_20260927.json` and
`metadata/historian_largest_family_floor_sensitivity_20260927.json`.

A combined largest-family comparison is queued behind both complete retry
audits. It retains the four combinations of MAFFT/FAMSA and minimum edge
length 1e-9/1e-7, requires identical ungapped inputs for all 622 extant proteins,
and checks candidate descendant sets. Sixteen candidate comparisons change
only one factor at a time; four involve the assumed root. The output retains
all diagnostic sequences but does not qualify an ancestral FASTA.
Implementation: `scripts/compare_historian_largest_family.py`, with matching
plan and launch metadata. One CPU, 4 GiB RAM, no swap, 0.1 GiB output and
0.01–0.5 hours after the audits are allowed. No new inference or GPU work
is launched by this comparison. The process is verified live and waiting.

## Complete four-alternative largest-family comparison

Both FAMSA retries and their independent audit completed successfully. Their
non-root branch-floor edit distances are 8, 14 and 19 at levels 0–2; the
assumed-root distance is 16. Combined with the completed MAFFT retries, all
16 one-factor candidate comparisons change sequence. Across the six non-root
alignment comparisons, edit distances range from 8 to 22; across the six
non-root floor comparisons, they range from 8 to 19. Four additional
comparisons concern the assumed root and remain marked separately.

All source and output hashes were checked, and a separate rolling dynamic
program reproduced all 16 edit distances from the retained sequences. Both
alignments contain the same 622 ungapped extant proteins; candidate descendants
and same-floor tree hashes are identical. The differences therefore concern
reconstruction sensitivity, not differences in the input protein sampling.
Approximate profile histories, selected traceback and fixed model settings
still limit interpretation. No alternative is qualified as a historical
ancestor by these diagnostic results. Original memory failures remain in the
capacity-grid record. Completed evidence is in
`metadata/historian_largest_family_comparison_completed_20260927.json` and
`metadata/historian_famsa_memory_retry_audit_completed_20260927.json`.


### Complete polytomy-resolution diagnostic within the continuing full grid

A new independent readback snapshot checked 292 successful original-grid
runs, retained four memory-limit outcomes, and left 28 jobs pending. The
snapshot is immutable; final all-job audit remains queued separately.
The four largest-family failures already have successful independently
checked higher-memory retries, which are intentionally preserved separately.

All 180 OG0000230 configurations are now independently checked: six alignment
inputs, two minimum-edge floors and all 15 rooted binary resolutions of the
four-way branch. Candidate sequences are identical across the 5,040 paired
node comparisons changing resolution alone and the 360 changing floor alone.
These totals include the explicitly assumed root and repeated comparisons of
shared inputs. They are not counts of independent evolutionary observations.
A family/dataset/factor/root summary is produced by
`scripts/summarize_historian_sensitivity_readback.py`; provenance is recorded in
`metadata/historian_capacity_sensitivity_partial_v3_20260927.json`.

Across the entire original-grid snapshot, all 5,624 available paired candidate
comparisons have zero edits. This must not be reported as general robustness:
the separate OG0000972 higher-memory retries show non-root differences of
8–19 edits across floors and 8–22 across input alignments. Neither result
establishes convergence or posterior concentration; Historian's approximate
profile and traceback limitations still apply. No final ancestral ensemble
has been qualified by this diagnostic.
