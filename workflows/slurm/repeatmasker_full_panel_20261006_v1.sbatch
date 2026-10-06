#!/usr/bin/env bash
#SBATCH --job-name=fungal-rm-full
#SBATCH --partition=gpu
#SBATCH --array=0-525%12
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/repeatmasker-full-panel-20261006-v1/%A_%a.out
#SBATCH --error=logs/repeatmasker-full-panel-20261006-v1/%A_%a.err
set -euo pipefail

cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/STRUCTURAL_GENOMICS
mkdir -p logs/repeatmasker-full-panel-20261006-v1
ENV=.cache/envs/fungal-genome-architecture-repeat/bin
export PATH="$PWD/$ENV:$PATH"
"$ENV/python" scripts/run_repeatmasker_taxon_v1.py \
  --config config/repeatmasker-full-panel-20261006-v1.json \
  --task-index "$SLURM_ARRAY_TASK_ID" \
  --library-builder scripts/build_repeatmasker_library_v1.py \
  --repeatmasker "$ENV/RepeatMasker"
