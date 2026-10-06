#!/usr/bin/env bash
#SBATCH --job-name=fungal-repeat-proximity
#SBATCH --partition=gpu
#SBATCH --array=0-525%12
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/repeat-proximity-full-panel-20261006-v1/%A_%a.out
#SBATCH --error=logs/repeat-proximity-full-panel-20261006-v1/%A_%a.err
set -euo pipefail

cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/STRUCTURAL_GENOMICS
mkdir -p logs/repeat-proximity-full-panel-20261006-v1
ENV=.cache/envs/fungal-genome-architecture-repeat/bin
export PATH="$PWD/$ENV:$PATH"
"$ENV/python" scripts/run_repeat_annotation_proximity_taxon_v1.py \
  --config config/repeat-proximity-full-panel-20261006-v1.json \
  --task-index "$SLURM_ARRAY_TASK_ID" \
  --parser-script scripts/parse_repeatmasker_taxon_v1.py \
  --linker-script scripts/link_gene_repeat_proximity_v1.py
