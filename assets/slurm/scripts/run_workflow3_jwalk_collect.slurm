#!/bin/bash
#SBATCH --job-name=workflow3_jwalk_collect
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=6:00:00
#SBATCH --output=assets/slurm/logs/workflow3_jwalk_collect-%j.log
#SBATCH --error=assets/slurm/logs/workflow3_jwalk_collect-%j.err

# Workflow 3 Step 4b: Collect Jwalk/XP data for all 1000 trajectories
# This job collects per-trajectory XP cross-link probability data into consolidated arrays
# Expected runtime: 3-4 hours

set -e

echo "=========================================="
echo "Workflow 3 Step 4b: Jwalk Collection"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo ""

# Activate environment
source ~/.bashrc
conda activate entdetect

# Run collection
python collect_jwalk_optimized.py

echo ""
echo "Completion time: $(date)"
echo "=========================================="
