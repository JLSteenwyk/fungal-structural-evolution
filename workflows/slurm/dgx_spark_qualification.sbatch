#!/bin/bash
#SBATCH --job-name=fungal-esmfold-check
#SBATCH --partition=spark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/qualification-%j.log
set -euo pipefail
cd /home/jlsteenwyk/projects/fungal-structural-evolution
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
exec esmfold-env/bin/python scripts/run_marker_predictions.py \
  --inputs data/qualification --checkpoint data/checkpoint \
  --output results/qualification-v1 --max-length 768 --limit 4
