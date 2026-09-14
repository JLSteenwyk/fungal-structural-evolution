# Temporary cooling for the September 14 runs

User requested cooler operation for the current work, with slower completion
acceptable. This is an execution exception, not a new project default.

The active jobs were launched directly, outside SLURM. Editing a batch script
would not change them. Setting GPU 0 to 200 W failed because administrator
permission was unavailable; the hardware limit remains 300 W.

The temporary controller restricts the recorded same-user project processes
and their descendants to 24 logical CPUs representing 24 distinct physical
cores, shared across the jobs. All their existing threads are restricted;
new threads and descendants are checked every two seconds. It also alternates
30 seconds of execution and 30 seconds of suspension for the exact live
ESMFold process (PID plus creation time). Already submitted GPU work may
briefly continue after suspension. This lowers average activity, not peak
power or memory use. Predictions, models and scientific settings are unchanged.

The controller expires after 48 hours, approximately September 16 at 00:56
Eastern. SIGTERM also resumes its suspended GPU process and restores recorded
CPU affinities. No startup service, persistent GPU cap or SLURM default was
installed. Current configurations are recorded in
`metadata/temporary_cooling_20260914_config.json` and launch identity in
`metadata/temporary_cooling_20260914_launch.json`. Do not reuse the recorded
process IDs for another run.

To end this particular adjustment early, verify the controller identity in
that launch record, then send SIGTERM to the controller (PID 70110 at launch).
Do not kill the ESMFold process. If the controller is killed abruptly, recovery
is available after checking it is no longer alive:

```bash
python scripts/temporary_run_resource_limit.py \
  --config metadata/temporary_cooling_20260914_config.json \
  --state results/temporary_cooling_20260914_state.json --restore
```

Validation: a synthetic process test verified restricted CPU affinity,
suspension, resumption and affinity restoration on SIGTERM. Initial live
readback checked 913 thread affinities with no mismatches. GPU samples during
an initial rest period fell from roughly 298 W / 81–82 C under load to
83–89 W / 69–74 C. These are transient samples, not steady-state thermal or
power guarantees. ESMFold output generation is checked separately. The former
41-hour ETA no longer applies; half-time GPU scheduling alone approximately
doubles remaining prediction time, with further CPU contention possible.

## Temporary full-compute interval, September 14

At the user's request, the original cooling controller was stopped cleanly
at about 09:10 Eastern. Full affinity restoration passed 896 thread checks;
ESMFold resumed uninterrupted execution with access to all 192 logical CPUs.
The prior controller PID 70110 is no longer active.

A user systemd timer, `fungal-reapply-cooling-20260914.timer`, is armed for
September 14 at approximately 12:39:56 Eastern, 3.5 hours after scheduling.
It will start the same resource controller with the same 24-core and 30-second
run/rest settings, using the recorded surviving project process trees and
their descendants. It retains the original September 16 00:56 Eastern expiry.
This one-off timer changes no future analysis defaults and does not require
this chat to remain active; a host reboot would remove the transient timer.

Schedule, exact commands and checks:
`metadata/temporary_full_compute_20260914_schedule.json`.
Reapplication configuration:
`metadata/temporary_cooling_20260914_reapply_config.json`.
After activation, the new controller's state and recovery information will be
in `results/temporary_cooling_20260914_reapplied_state.json`. Use these new
paths, rather than the original state, for subsequent recovery. The systemd
service is `fungal-reapply-cooling-20260914.service`; stopping that service
normally sends SIGTERM and restores the recorded settings.
