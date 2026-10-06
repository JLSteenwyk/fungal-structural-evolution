#!/usr/bin/env bash
#SBATCH --job-name=fungal-gene-anchors
#SBATCH --partition=gpu
#SBATCH --array=0-525%12
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/gene-product-anchors-full-panel-20261006-v1/%A_%a.out
#SBATCH --error=logs/gene-product-anchors-full-panel-20261006-v1/%A_%a.err
set -euo pipefail

cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/STRUCTURAL_GENOMICS
mkdir -p logs/gene-product-anchors-full-panel-20261006-v1
/home/bizon/anaconda3/bin/python scripts/run_gene_product_anchor_taxon_v1.py \
  --config config/gene-product-anchors-full-panel-20261006-v1.json \
  --task-index "$SLURM_ARRAY_TASK_ID"
