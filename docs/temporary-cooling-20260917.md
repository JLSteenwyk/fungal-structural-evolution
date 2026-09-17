# Temporary reduced-heat execution, September 17

At the user's request, the current project jobs have returned to the earlier
30-second run / 30-second rest pattern, now applied to both ESMFold workers.
Recorded project process trees share the previous 24-core CPU set. This is a
one-off execution adjustment, not a new project or SLURM default. Scientific
prediction settings and GPU hardware power limits remain unchanged; the
rest periods reduce average activity rather than peak hardware power.

Two instances of the existing tested `temporary_run_resource_limit.py` control
nonoverlapping process trees. GPU 0's controller also limits the current CPU
analyses; GPU 1's controller owns its prediction wrapper and descendants.
The launch record captures each PID and creation time. The controller saves
original thread affinities before modification and resumes its GPU worker
and restores affinities on normal expiry or SIGTERM.

The deadline is September 17 at 18:00 America/Los_Angeles (PDT), equivalent
to September 18 at 01:00 UTC. The user said PST; Pacific local time is the
stated assumption pending any clarification. Literal 18:00 PST would instead
be 19:00 PDT. `fungal-end-cooling-20260917.timer` invokes a normal stop of both
controllers at the deadline. Each controller also has its own duration bound.
The timer and services are transient and do not persist across a host reboot.

Configuration and process identities are in
`metadata/temporary_cooling_20260917_launch.json` and the two corresponding
GPU config files. The live restoration states are
`results/temporary_cooling_20260917_gpu0_state.json` and
`results/temporary_cooling_20260917_gpu1_state.json`. Stop the controller
services normally to end this adjustment early; do not terminate the ESMFold
workers. The documented `--restore` mode of the original controller remains
available for recovery after verifying a controller is no longer alive.

Initial readback checked 90 live thread affinities with no mismatches. Run/rest
cycle verification and the timer readback are recorded separately in
`metadata/temporary_cooling_20260917_verification.json`. Prior GPU completion
estimates must be adjusted for this temporary reduction in active time.
