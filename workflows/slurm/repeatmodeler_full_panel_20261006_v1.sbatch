#!/usr/bin/env bash
#SBATCH --job-name=fungal-rm2-full
#SBATCH --partition=gpu
#SBATCH --array=0-525%12
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/repeatmodeler-full-panel-20261006-v1/%A_%a.out
#SBATCH --error=logs/repeatmodeler-full-panel-20261006-v1/%A_%a.err
set -euo pipefail

cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/STRUCTURAL_GENOMICS
mkdir -p logs/repeatmodeler-full-panel-20261006-v1
ENV=.cache/envs/fungal-genome-architecture-repeat/bin
export PATH="$PWD/$ENV:$PATH"
"$ENV/python" scripts/run_repeatmodeler_taxon_v1.py \
  --config config/repeatmodeler-full-panel-20261006-v1.json \
  --task-index "$SLURM_ARRAY_TASK_ID" \
  --build-database "$ENV/BuildDatabase" \
  --repeatmodeler "$ENV/RepeatModeler"
