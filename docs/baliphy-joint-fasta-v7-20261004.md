# Joint ancestral alignment serialization, October 4

The current candidate writes the same joint ancestral-state record using
sequence rows directly, avoiding conversion of the entire multi-megabyte FASTA
alignment into a Haskell character list and then back into text. The ancestral
sequence calculation, sampled alignment, substitution/indel models, prior,
tree and native resource limits remain unchanged. This is a candidate software
correction; the full failure grid and ancestral posterior are not accepted.

A diagnostic wrapper first marked the start/end of each logger's context and
write actions. All 405 full model programs reverse exactly after removing the
markers. Three existing native prior controls preserve six output files byte
for byte, with 63 independently read scalar rows and 1,008 phase markers.
[Original software wait and complete journal proof](../metadata/baliphy_logger_phase_trace_software_transport_20261004_v1.json).

The separate full-input phase diagnostic captures SIGSEGV during iteration-zero
joint-state writing. The fault is eight bytes below the fully grown 8 MiB stack,
at a native evaluator push instruction. No joint frame or scalar row is saved.
Original wait 64683 and complete wrapper/manager records close; debugger exit
zero records a captured native fault, not a completed sampler.
[Fault evidence](../metadata/baliphy_logger_phase_failure_review_20261004_v1.json).

The original **8 MiB uninstrumented** debugger output also has no saved scalar
rows or joint frames. The distinct **64 MiB** uninstrumented comparison saves
iterations zero through nine before native SIGKILL. The phase review's V1
statement about later uninstrumented iterations must be read with that
distinction. Phase markers can change evaluation context; neither tiny-control
noninterference nor this one signal stop attributes every original failure to
a particular field or proves a repair.

The V7 candidate obtains the existing ancestral character data once, preserving
the original FASTA writer, and constructs each joint header and sequence line
from the same labels, alphabet, ambiguity database and native sequence vectors.
It replaces only the `Text.unpack` / `lines` / `Text.pack` roundtrip. All 405
model sources reverse exactly; the broad, centered and package native controls
preserve TSV, corrected scalar JSON, column map, tree, ancestral FASTA and joint
JSON byte for byte. All nine native joint frames and 2,709 scalar comparisons
pass. Installed software, priors, tolerances and original output bytes are
unchanged.
[Candidate controls and original transport](../metadata/baliphy_joint_fasta_v7_software_transport_20261004_v1.json).

A fresh standalone comparison uses the original failed family's **622 tips and
1,243 tree nodes**, input, broad prior, seed 1322110685 and 20 iterations. Its
native limits remain **48 GiB address space, 5,674 CPU seconds, 2 GiB per file,
8 MiB soft/unlimited hard stack**. No debugger or phase markers are used.
Wrapper/native share two CPUs, 64 GiB memory, no swap and one BLAS thread. The
prospective 2–15 minute estimate is uncalibrated, with the original hard caps
retained; it is not a completion ETA. The local output allowance is 4 GiB and
no paid infrastructure or GPU is used.

At 14:13 UTC, original wait **58381** is live, with the native command and
actual caps verified. A closed observation snapshot contains six scalar rows
and the first 7,519,894-byte joint frame. Those nonempty bytes match the earlier
64 MiB run exactly. Independent readers check 258 scalar values, 622 tips,
621 ancestors and 320,014 ancestral residue/category pairs. Empty files from
the original 8 MiB diagnostic would be vacuous prefix evidence and are excluded
from this compatibility claim.
[Closed prefix snapshot, source hashes and live identity](../metadata/baliphy_joint_fasta_v7_prefix_checkpoint_20261004_v1.json).

The comparison root is
`results/ancestral/baliphy-joint-fasta-v7-full-input-comparison-20261004-v1/`;
the closed prefix lives under
`data/software_audits/baliphy-joint-fasta-v7-prefix-observation-20261004-v1/`.
Large/raw outputs stay outside Git and are bound by the receipts. Key scripts
are `baliphy_joint_fasta_lines_v7.py`, `check_baliphy_joint_fasta_lines_v7.py`,
`prepare_baliphy_joint_fasta_v7_full_comparison_v1.py`,
`run_baliphy_joint_fasta_v7_full_comparison_v1.py` and
`record_baliphy_joint_fasta_v7_prefix_checkpoint_v1.py`.
[Frozen execution plan](../metadata/baliphy_joint_fasta_v7_full_comparison_plan_20261004_v1.json).

Continue the same original wait. A completed native horizon must receive full
independent output readback and original transport closure before comparison
across all 24 failed roles. Adequate posterior sampling and all eight scientific
aims remain separate unfinished requirements. No original failure is retried,
replaced or interpreted as successful, and no family or fungal lineage is
dropped to obtain a passing run.

## Completed full-input comparison and all 24 original failed roles

The original standalone native sampler has now completed all 20 iterations,
with exit zero and 21 scalar rows/903 mapped values/three joint frames. Runtime
was 1,429.61 seconds at the original caps. The wrapper subsequently exited one
while looking for the column map in the earlier empty 8 MiB diagnostic. This
reporting error, source, output namespace and original wait remain preserved;
a fresh V2 read-only review uses the nonempty 64 MiB prefix and does not rerun
the native attempt. All six older nonempty files match the new output prefix.
[Exact original native and reporting outcomes](../metadata/baliphy_joint_fasta_v7_full_comparison_review_20261004_v2.json).

A separate independent reader checks every native topology/branch, scalar,
FASTA, joint-state frame and projection over 622 tips and 621 ancestors.
All 903 scalar values, 956,096 ancestral residue/category pairs, 854,976 tip
pairs and twelve projected candidate frames pass. Original reader wait 42976,
whole wrapper payloads and invocation-bound manager records close successfully.
This proves the output integrity of one controlled full-input horizon, not
adequate posterior sampling or a repair across every original failure.
[Independent full readback](../metadata/baliphy_joint_fasta_v7_full_comparison_readback_20261004_v1.json).

The new all-failure comparison covers **24 original roles, six prior quartets,
two effective full inputs and all three priors**. Software qualification
checks each actual original recipe and rejects twenty altered contracts.
Original seeds, inputs, priors, CPU/wall/file/address-space and 8 MiB stack
limits are retained. The only inference-program change is the qualified V7
serialization transformation. Old failed attempts are never overwritten.
Four workers hold 48 GiB address-space leases under a 200 GiB/no-swap cgroup;
leases are reservations, not measured peak usage. The one completed 24-minute
case supports only an uncalibrated 2–12 hour planning range for this grid.
[Full plan](../metadata/baliphy_joint_fasta_v7_failure_grid_plan_20261004_v1.json).

At 15:35 UTC, four broad-prior roles from the first effective input have native
zero exits. Three pass finite-output checks; chain three retains an explicit
nonfinite/literal-null review and has no accepted array. These are producer
outcomes, pending whole-grid independent readback and original transport.
Four further workers remain live. Original producer wait 72605 is unchanged;
reader wait 84274 is queued with two CPUs/32 GiB/no swap behind original
producer invocation success. The reader reconstructs every outcome with
export creation disabled. Neither native zero exits nor review-state tagging
qualifies an ancestral posterior. All eight scientific aims remain open.
[Actual runtime](../metadata/baliphy_joint_fasta_v7_failure_grid_checkpoint_20261004_goal_1535.json).

At 15:47 UTC, the original timing stage reaches fifteen producer checkpoints
with sixteen actual workers and no failure files. The all 24 ancestral
comparison has seven completed producer dispositions: six finite-output
checks and one retained special-value review, all native exits zero.
Original producer wait 72605 and reader wait 84274 remain live; the reader has
no reconstruction child before producer closure. These partial outcomes do
not qualify the entire failure grid or an ancestral posterior.
