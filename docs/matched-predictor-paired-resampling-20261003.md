# Paired predictor branch uncertainty

The full paired resampling control is running. It extends the
[completed matched point fits](matched-predictor-branch-controls-20261003.md)
with uncertainty from sampled alignment columns. It does not replace the
approximately 500-fungus/25-outgroup study or establish accepted biological
effects.

## Scope and conditioning inputs

Use every one of the **133 ready inputs** from the closed 125-marker/70-view
grid. All 8,750 comparisons remain accounted upstream, including 3,780
insufficient comparisons. The usable predictor overlap remains 71 markers
and 21 fungal taxa, heavily selected toward Suillus. Each fit contains 4–16
taxa. It contains no outgroup and is not representative of fungal diversity.

For each input, run **200 independent-site draws and 200 circular-block
draws**, each the same length as the original alignment. Blocks contain ten
adjacent retained columns and wrap at the alignment boundary; trim the last
block to the original length. These are retained-column positions, which
can have gaps in original coordinates. Ten columns is a declared sensitivity
setting, not an established biological correlation length.

Apply each draw to all taxa and to AA, AlphaFold states and ESMFold states
identically. Keep the original eligibility, unknown masks, predicted states,
models and fixed topology. Do not drop taxa or reject and redraw a replicate
because its sampled observations become poor. Failed/invalid fits retain
their original identities and empty numerical slots.

Seventy-one groups with identical marker, taxa, source columns and all three
alignment hashes share draws across alternative topologies. This preserves
paired topology comparisons; alternative trees do not become independent
biological replicates. Seeds use SHA256-derived 128-bit integers and an
explicit NumPy PCG64 generator. Case identities include the original input,
mode, replicate and source-index hash. Predictor counterparts also share
the model-specific native optimization seed.

Each draw has seven native roles: AA LG+F+G4, plus both predictors under
AF+G4, AF+F+G4 and LLM+G4. There are **53,200 cases, 372,400 native fits and
6,146,000 branch-value slots**. IQTree preserves the supplied topology and
identical sequences. [Official command reference](https://iqtree.github.io/doc/Command-Reference).

## Qualification and audit evidence

The software gate independently reproduced all 53,200 draw identities and
site/block index arrays, then checked complete serialized draw accounting.
Fourteen successful synthetic native fits cover both modes and all seven
roles. Seven deliberately CPU-zero native fits check retention of real
failed processes. Artificial failed/review arrays preserve NaNs; sixteen
malformed grid, index, shape/dtype/value and archive contracts are rejected.
These are software checks within the full experiment, not a biological pilot.
[Software proof](../metadata/matched_predictor_resampling_software_validation_20261003_v1.json).

The bounded software service exited zero and recorded 22.61 child CPU seconds
and 305,098,752 bytes peak child RSS. Original invocation journals and actual
wait return codes are retained. Manager resource summaries are separate from
native process RSS. [Execution](../metadata/matched_predictor_resampling_software_execution_20261003_v1.json)
and [original wait evidence](../metadata/matched_predictor_resampling_software_transport_20261003_v1.json).

Every production case retains a lossless gzip tar archive of inputs, native
configuration, PID/create/command records, logs, checkpoints, reports, trees,
failed evidence and role dispositions. Every archived file is byte-hashed and
reread before deleting only that case's newly declared scratch directory.
Compact NPZ files retain sampled indices, seven status codes and all branch
values. Valid roles contain finite nonnegative values; failed/review roles
contain NaNs. Archive extraction rejects foreign paths, links, duplicate
members and missing/altered files.

The independent full reader reconstructs all resampled FASTA/configuration
files, replays original native command identities and checks every archive
member. DendroPy independently decodes native tree/report branch graphs;
all NPZ keys, shapes, dtypes and values are checked against archived roles.
Full producer/reader/source/artifact closure requires both exact original
invocation journals. Numerical output integrity does not establish model
adequacy or accepted evolutionary effects.

At **15:50 UTC on October 3**, 575 immutable case receipts were available.
A read-only early check covered six completed cases across both modes and
two original inputs: all 42 native roles passed independent archive/tree/
report/array readback. One live native worker's actual 2-GiB address-space,
300-CPU-second and 8-MiB file limits were verified. This is explicitly partial;
the complete reader remains queued behind the original full producer.
[Early production readback](../metadata/matched_predictor_resampling_partial_readback_20261003_v1.json).

## Resources and reproducibility

The original producer uses eight CPUs, eight workers, 24 GiB RAM and no
cgroup swap. Each case executes its seven roles sequentially, so at most
eight native roles run simultaneously. Each native role uses one thread,
2 GiB address space, 300 CPU seconds, 600 wall seconds and an 8-MiB file limit.
Reader and closure use two CPUs/16 GiB/no swap. No GPU prediction or new
cost is incurred; existing ancestral and covariance jobs stay unchanged.
[Frozen plan and resource estimate](../metadata/matched_predictor_resampling_plan_20261003_v1.json)
and [original controllers](../metadata/matched_predictor_resampling_launches_20261003_v1.json).

The closed point-fit directory contained 11,177 files/36,023,503 bytes.
Scaling its bytes by 400 gives approximately 13.4 GiB of raw output. Native
worker-wall scaling gives 5.67 hours at eight workers, **excluding controller,
input, archiving, reader overhead and changed optimization behavior**. This
is not a finish ETA. The output planning allowance is 128 GiB; a 256-GiB
free-disk guard is checked initially and at completed cases. These planning
figures are separate from the enforced per-process/file/cgroup limits.

Sources are `matched_predictor_resampling.py`,
`run_matched_predictor_resampling.py`,
`readback_matched_predictor_resampling.py` and
`prepare_and_launch_matched_predictor_resampling.py`. All source, input,
native binary and model hashes are frozen in the plan. The
[environment](../environments/matched-predictor-branch-control-20261003.yml)
records dependencies. `record_matched_predictor_resampling_checkpoint.py`
requires a new output filename and checks the exact original controllers.
Reproduction needs fresh explicit plan/output identities and updated pins;
never overwrite completed attempts or run an expired GPU window.

## Interpretation still required

After complete numerical/provenance closure, summarize paired branch changes
and compare site/block distributions. An unresolved replicate prevents
accepting an interval based only on successful fits. Two hundred draws give
modest tail resolution; any confidence-level claim also needs Monte Carlo
and coverage assessment. No intervals or accepted acceleration are reported
from partial production output.

Predicted structures, topology, orthology and confidence-selection remain
conditioning inputs. This does not resample biological structure-prediction
error, alignment uncertainty or full phylogenetic uncertainty. Branch lengths
remain state substitutions/site, not physical displacement or rates/year.
Boundary calibration, direct-coordinate/experimental controls and broad
fungal replication remain required. All eight biological aims are incomplete.
