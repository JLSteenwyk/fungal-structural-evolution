# Full-input latent-alpha logger comparison

The complete **24-role V10 comparison** launched at **23:42 UTC on October 5**.
Its original tool session is **59687**, and the exact wrapper/native/cgroup
checkpoint confirms four live native samplers. No role has completed at the
23:43 UTC checkpoint. Full native completion, paired scientific-file checks
and independent latent-state reconstruction remain required.
[Frozen plan](../metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json),
[original launch](../metadata/baliphy_log_alpha_v10_full_grid_launch_20261005_v1.json)
and [exact runtime observation](../metadata/baliphy_log_alpha_v10_full_grid_checkpoint_20261005_goal_launch_v1.json).

This covers both original 622-protein-tip effective inputs of OG0000972,
all three priors and four chains per prior: six complete quartets. Protein
tips are not independent species. The original 24 seeds, inputs, priors,
20 iterations, 48 GiB native address-space limits, 8 MiB stack and individual
CPU/wall/file limits are preserved. Only the generated logger program and
the fresh output namespace change. Original V6 failures and V7 special-value
reviews remain retained; this does not replace or continue an old chain.

This is the complete affected native comparison matrix within the project's
full workflow. It does not establish posterior convergence or shrink the
501-fungal-entry plus 25-outgroup study to one family. All eight evolutionary
aims, accepted frameworks and adequate ancestral uncertainty remain incomplete.

## Source and software prerequisites

The [V10 logger](baliphy-latent-log-alpha-v10-20261005.md) already passes all
405 generated-source reversibility checks, three paired native prior controls,
18 byte-identical scientific-file comparisons and separate reconstruction of
120 diagnostic rows. Its successful software and original execution records
remain immutable. Earlier V8/V9 failures remain separate.

The new full job constructor checks every role against the complete closed V7
matrix and actual original V6 failed-role source records. All six generated
V10 programs reverse exactly to their original V7 source. Twenty deliberate
changes to roles, inputs, seeds, priors, source identities, native arguments,
memory reservations and limits are rejected. The original software session
3658 exits zero, with complete native/journal/source closure.
[Construction checks](../metadata/baliphy_log_alpha_v10_full_jobs_software_validation_20261005_v1.json).

Full-reader integration checks reuse all 63 already closed native V10 prior
control rows. Four deliberately changed/incomplete states remain excluded:
changed density, a 20-row partial trace, a native failure despite a complete
trace, and missing diagnostic output. Original software session 83547 exits
zero with complete native/journal/source closure. No new model run or corpus
pilot was used for these software checks.
[Reader integration checks](../metadata/baliphy_log_alpha_v10_full_reader_software_validation_20261005_v1.json).

Source preparation completes under its actual original tool session 64935.
The immutable plan binds **20,210 sources**, both software closures, all
original native pairs, the original model/input/configuration chain, the new
programs and the complete producer/reader code. Preparation is separate from
the actual native launch.

## Resources and observed baseline

The producer uses four CPU workers, 200 GiB RAM, no swap and one BLAS thread.
Four 48 GiB native leases reserve 192 GiB; 8 GiB remain for the controller.
The controller keeps an 8 GiB soft/192 GiB hard address-space limit so each
native child can retain its original 48 GiB limit. The stage wall cap is
24 hours, with individual native limits preserved and no automatic retry.
The 128 GiB output allowance and 256 GiB minimum free-disk reserve are
prelaunch estimates, not current output size or measured native memory.
[Producer resources](../metadata/baliphy_log_alpha_v10_full_grid_resources_20261005_v1.json).

The previous complete V7 matrix consumed 10.443 aggregate worker-hours;
individual native runs took 1,179–1,947 seconds. Dividing the worker total
by four gives a 2.611-hour lower bound that omits scheduling, validation and
V10 overhead. The V10 planning range is **2–12 hours**, explicitly uncalibrated
for this new version. Native limits and planning ranges are not ETAs or
posterior sampling requirements. No GPU, paid provisioning, change to another
project or installed-library edit is involved.

## Every outcome and the independent reader

The producer retains all 24 native/scalar/joint dispositions, including native
failures and explicit special values. It compares six original scientific
files per role: scalar TSV/JSON, column map, runtime tree, ancestral FASTAs
and joint site-property samples. All **144 planned file comparisons** retain
their original/current paths and hashes; missing or changed files produce
review outcomes rather than dropped roles. The new diagnostic trace is kept
separately in each native attempt.

The separate reader is prepared, not launched or queued. Actual original full
producer API/native zero, complete invocation-journal matching and all source
bindings must close before its launch. It independently reconstructs every
file comparison with a separate byte reader, then applies the qualified strict
JSON/90-digit Decimal diagnostic reader to every available latent row. It
checks the saved context, prior parameters/density, derived alpha, rate state,
quality tags and complete 21-row trace. Partial/native-failed/malformed traces
remain explicit. Per-row numeric errors are retained without silently dropping
the rest of the trace. Existing native/scalar/joint inspectors and job contracts
are shared dependencies and are stated as such.

The reader uses two CPUs, 32 GiB RAM, no swap, 24 GiB hard/8 GiB soft
controller address space and one BLAS thread. Its 0.1–3 hour estimate is
uncalibrated; three-hour wall/two-hour CPU safety caps are not ETAs. It does
not launch native samplers or admit old reviewed arrays.
[Reader resources](../metadata/baliphy_log_alpha_v10_full_grid_readback_resources_20261005_v1.json)
and [full reader](../scripts/readback_baliphy_log_alpha_v10_full_grid_v1.py).

Even a successful full-input logger comparison will not make the special-value
arrays suitable for posterior inference. New latent captures can establish
the mechanism for newly observed states; they cannot recover unsaved historical
values. Adequate inference still requires estimated full sampling resources,
convergence, uncertainty propagation and calibrated biological interpretation.

## Reproduction

Large outputs remain outside Git at
`results/ancestral/baliphy-log-alpha-v10-all24-full-input-comparisons-20261005-v1/`.
Versioned artifacts include scripts, resource plans, original tool payloads,
software execution logs, launch identity and source hashes. Live native logs
are not committed. Complete raw-input public release remains outstanding.

Use fresh output/receipt namespaces with the pinned native installation and
original inputs. The full constructor and reader software checks must close
before source preparation; the exact original producer completion must close
before readback. Do not reuse an incomplete attempt or silently restart a
missing process handle. The producer command for this frozen plan is:

```bash
python scripts/run_baliphy_log_alpha_v10_full_grid_v1.py \
  --plan metadata/baliphy_log_alpha_v10_full_grid_plan_20261005_v1.json \
  --receipt metadata/baliphy_log_alpha_v10_full_grid_20261005_v1.json
```

The command must run through the recorded bounded wrapper and exact cgroup
resources. A clone alone does not include the large raw sources.
