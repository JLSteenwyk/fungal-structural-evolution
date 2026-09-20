# Prediction pause and overnight window, September 18–19

**Latest update, September 20:** prediction resumed at 06:53 PDT for six
hours, with automatic pause at 12:53:17 PDT today. See
`metadata/prediction_window_20260920_schedule.json`. Earlier windows below
are historical. The same process and immutable manifest remain in use.

**Updated schedule:** the user superseded the original window. Structure
prediction remains paused until **03:00 PDT September 19**, runs for ten
hours, then pauses at **13:00 PDT September 19**. The old two timers were
stopped and replaced by `fungal-prediction-resume-revised-20260919.timer` and
`fungal-prediction-pause-revised-20260919.timer`. Both are verified active.
The immutable original manifest/initial pause state are reused for resume;
the final pause writes `results/prediction_window_20260919_revised_final_pause.json`.
See `metadata/prediction_window_20260919_revised_schedule.json`.

The earlier schedule below is retained as history and is no longer active.

The user requested a complete pause of protein structure prediction until
18:00 Pacific local time on September 18, followed by 12 hours of execution
and another pause at 06:00 Pacific local time on September 19. PDT is used,
consistent with the user's prior explicit timezone clarification.

The only active prediction worker, GPU 1 ESMFold PID 267239 (creation time
1789573381.2), was suspended in place with SIGSTOP. Completed models and live
inference state are preserved. GPU 0's previous batch had already finished;
no new GPU worker was started. CPU phylogeny/domain analyses are unaffected.

The existing identity-checking `pause_resume_project.py` script handles both
scheduled actions. The resume timer is `fungal-prediction-resume-20260918.timer`
(September 19 01:00 UTC); the subsequent pause timer is
`fungal-prediction-pause-20260919.timer` (September 19 13:00 UTC). The new pause
uses a separate state file. Both timers were verified active. These are
transient user systemd timers, so a host reboot requires explicit schedule
review and prediction recovery. The chat need not remain active.

The manifest, schedule and verification are recorded under
`metadata/prediction_window_20260918_*`. The initial live state is
`results/prediction_window_20260918_initial_pause.json`; the next pause state
will be `results/prediction_window_20260919_final_pause.json`. Keep the manifest
and controller script unchanged while scheduled. Any later queue rebalance
must preserve this requested window and update process identities/timers;
these timers intentionally do not signal unrelated replacement PIDs.

Inference elapsed times spanning the suspension include paused time and must
not be treated as full-speed benchmarks. This schedule resumes the existing
single-GPU queue; it does not implement the discussed two-GPU rebalance.
