#!/usr/bin/env bash
# Submit only after Slurm is reachable and the original array 23002 is terminal.
#SBATCH --job-name=fungal-rm2-recovery
#SBATCH --partition=gpu
#SBATCH --array=11,37,102,143,154,156,188,189,202,242,249,274,276,290,364,378,409,522,523,524,525%12
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/repeatmodeler-recovery-20261021-v1/%A_%a.out
#SBATCH --error=logs/repeatmodeler-recovery-20261021-v1/%A_%a.err
set -euo pipefail

cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/STRUCTURAL_GENOMICS
mkdir -p logs/repeatmodeler-recovery-20261021-v1
ENV=.cache/envs/fungal-genome-architecture-repeat/bin
export PATH="$PWD/$ENV:$PATH"
export FAMDB_DIR="$PWD/.cache/envs/fungal-genome-architecture-repeat/share/famdb-3.0.0"
"$ENV/python" scripts/run_repeatmodeler_taxon_v1.py \
  --config config/repeatmodeler-full-panel-20261006-v1.json \
  --task-index "$SLURM_ARRAY_TASK_ID" \
  --build-database "$ENV/BuildDatabase" \
  --repeatmodeler "$ENV/RepeatModeler" \
  --repeatclassifier "$ENV/RepeatClassifier" \
  --famdb-dir "$FAMDB_DIR"
