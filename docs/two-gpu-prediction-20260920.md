# Two-GPU structure prediction, September 20

**Latest schedule:** both workers resumed at 14:38:37 PDT September 20
for one hour, pausing at 15:38:37 PDT. They resume again at **23:30 PDT
September 20** and pause at **06:00 PDT September 21**. The three active
calendar timers are `fungal-two-gpu-pause-1h-20260920.timer`,
`fungal-two-gpu-resume-night-20260920.timer`, and
`fungal-two-gpu-pause-morning-20260921.timer`. The first pause state is
`results/prediction_two_gpu_20260920_one_hour_pause_state.json`; the final
overnight pause state is `results/prediction_two_gpu_20260921_overnight_pause_state.json`.
The same two workers, disjoint inputs and immutable pause manifest are used.
The earlier deadline below is retained as history.

The user authorized use of GPU 0 until the existing 12:53 PM PDT pause.
The original GPU 1 worker was asked to stop cleanly after its current model.
It exited with a recorded intentional interruption, no OOM, 2,318 completed
models and 2,045 unpredicted sequences. All completed model artifact hashes
were checked and the existing outputs remain unchanged.

`scripts/split_remaining_prediction_queue.py` created disjoint fresh inputs
under `data/prediction_inputs/markers-769-1024-two-gpu-20260920/`. GPU 0 has
1,022 sequences and GPU 1 has 1,023. Cubic sequence-length costs balance the
queues. An independent Bio.SeqIO readback verified every sequence against the
original FASTA, every marker link against the original table, disjointness,
and exact coverage of all 4,363 candidates by the completed set plus both
new queues. The six previously reusable models remain tracked separately
in the original length-band provenance.

Both workers use the unchanged ESMFold checkpoint and inference script,
full sequences, maximum length 1,024, chunk size 64, four CPU threads and
normal GPU power limits. Each service has a 64 GiB host-memory limit with
swap disabled. Both identical 48 GiB GPUs were checked idle before launch.
Existing long-protein timing observations give a planning range of 25–30
hours for the full remaining split at uninterrupted dual-GPU speed, but this
execution pauses at the user's much earlier deadline. No new charges.

The new outputs are
`results/predictions/esmfold-markers-769-1024-split-gpu0-20260920` and
`results/predictions/esmfold-markers-769-1024-split-gpu1-20260920`.
The two services and exact worker identities are recorded in
`metadata/prediction_two_gpu_20260920_launch.json`. The timer
`fungal-split-prediction-pause-20260920.timer` pauses both workers at
**12:53:17 PDT September 20**, writing
`results/prediction_two_gpu_20260920_pause_state.json`. The former single-worker
pause timer was canceled only after this replacement was verified active.

The original single-GPU wrapper reports its deliberately interrupted tier
as incomplete; it is superseded, not to be restarted. Future progress and
structure auditing must combine the 2,318 original completed models with the
two disjoint new output directories while preserving their distinct configs.
Use the new pause manifest and state for any authorized future resume; old
PID 267239 is no longer the predictor. Full output audits, conversion and
cross-source aggregation remain separate downstream work. Do not edit live
input receipts, configurations or the inference script.
