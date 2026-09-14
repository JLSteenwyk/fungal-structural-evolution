# DGX Spark resource assessment

The connected DGX Spark is accessible through the existing `spark` SSH alias
and direct 10 Gb/s Ethernet link. Observed September 14: GB10 GPU idle, 20 ARM
CPU cores, approximately 116 GiB available memory and 3.2 TiB free disk space.
Its SLURM `spark` partition is operational. A two-CPU, 8-GiB, one-GPU smoke job
completed a checked CUDA matrix calculation. No changes to the workstation's
scheduled cooling controls were made.

Spark has an ARM64 processor and shared CPU/GPU memory; the workstation's x86
binaries cannot simply be copied over. NVIDIA documents the architecture and
software requirements in its [porting guide](https://docs.nvidia.com/dgx/dgx-spark-porting-guide/overview.html)
and [software stack](https://docs.nvidia.com/dgx/dgx-spark-porting-guide/porting/software-requirements.html).

An isolated project virtual environment was created under
`/home/jlsteenwyk/projects/fungal-structural-evolution/esmfold-env` on Spark.
It inherits PyTorch 2.10.0+cu128 from the existing `glm_safety` environment;
Transformers 4.57.6 and dependencies were installed only in the project
virtual environment. The base environment was not modified. Exact relevant
versions are recorded in metadata/dgx_spark_environment_versions.json. This
is a documented environment with an external base dependency, not a fully
self-contained or immutable environment.

The existing PyTorch installation warns that its listed maximum CUDA
capability is 12.0 while GB10 reports 12.1. A successful simple GPU operation
alone does not establish full ESMFold compatibility. Qualification job 20912
therefore runs four already-predicted proteins of lengths 128, 256, 500 and
700, selected only by length, using the unchanged runner and pinned model
weights. It requests four CPUs, 64 GiB, one GPU and 30 minutes on Spark.
Source/checkpoint hashes are checked by the runner; outputs are stored
separately under `results/qualification-v1` on Spark. This is a hardware
compatibility check within the full project, not a biological pilot.

Job script: workflows/slurm/dgx_spark_qualification.sbatch. Resource estimates,
input selection and inventory are in metadata/dgx_spark_qualification_* and
metadata/dgx_spark_resource_inventory.json. The runner and its checksum helper
were copied unchanged; the 8.44-GB weights were transferred with symlinks
resolved. No production queue was split or moved during this assessment.

Appropriate next uses are disjoint prediction batches and ARM-compatible CPU
analysis tasks. The existing Spark HMMER installation reports version 3.4;
IQ-TREE and other workstation binaries need native builds and version checks.
Large species-tree fits already using hundreds of GiB should remain on the
workstation. Any production allocation must record a disjoint input manifest,
reserve Spark resources through SLURM, and audit returned files before they
enter the structural atlas. Hardware/software variants remain explicit in
prediction provenance and sensitivity analyses.

## Completed qualification and disposition

SLURM job 20912 completed in 5 minutes 6 seconds, exit 0:0. All four predictions
completed without OOM or interruption. Source artifacts, saved configurations,
complete sequences, finite coordinates and confidence arrays were read back.
Independent BioPython rigid-superposition calculations agree with the primary
RMSD calculation within 9.72e-16 Angstrom.

| Length | Workstation inference seconds | Spark inference seconds | Whole-chain CA RMSD, Angstrom |
| --- | ---: | ---: | ---: |
| 128 | 0.92 | 4.48 | 0.0833 |
| 256 | 3.26 | 14.99 | 0.3063 |
| 500 | 18.89 | 77.17 | 0.0048 |
| 700 | 42.37 | 175.54 | 6.4171 |

These four checks took approximately 4.1–4.9 times longer per prediction on
Spark under the tested configurations. This is not an exhaustive or controlled
hardware benchmark: platform, PyTorch version and available CPU resources
differ. At two sampled instants Spark reported approximately 55–57 W GPU
power while computing; these are not wall-power or energy measurements.

The 700-residue comparison remains materially different after restricting
both structures to CA pLDDT >=70: 376 positions yield RMSD 1.6456 Angstrom.
Whole-chain mean confidence differs little, illustrating why it is insufficient
as an equivalence check. The 128-residue control has no jointly >=70 positions
and should not be presented as a high-confidence match. These results do not
identify whether the discrepancy originates from software, platform arithmetic,
structural flexibility or another source, or establish which prediction is
more accurate. They do not calibrate a noise floor for evolutionary analyses.

Spark is available for additional computation, but its prediction outputs are
not approved for pooling with the workstation atlas. Resolve the longer-chain
platform discrepancy before production GPU assignment, or explicitly treat
Spark as a separate prediction variant with appropriate validation. Native CPU
workloads remain a separate route for added capacity. No production jobs were
launched; Spark is idle again after qualification.

Full comparison and confidence-mask readback:
metadata/dgx_spark_qualification_comparison.json and
metadata/dgx_spark_qualification_independent_readback.json. Scripts:
scripts/compare_spark_qualification.py and scripts/readback_spark_qualification.py.
The checkpoint and returned prediction files remain outside Git; source hashes,
versions and execution settings are versioned. Historical SLURM accounting
(`sacct`) failed because its database connection was unavailable; terminal
`scontrol` job state and complete project receipts were retained instead.
