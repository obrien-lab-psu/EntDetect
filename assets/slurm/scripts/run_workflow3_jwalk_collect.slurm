#!/bin/bash
#SBATCH --job-name=workflow3_jwalk_collect
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=5:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --output=assets/slurm/logs/workflow3_jwalk_collect-%j.log
#SBATCH --error=assets/slurm/logs/workflow3_jwalk_collect-%j.err

# Workflow 3: Jwalk/XP Collection SLURM Script
# 
# Collects per-trajectory Jwalk (cross-link probability) outputs into a 
# consolidated Jwalk.npy array.
#
# Expected runtime: 3-4 hours
# Memory: 32 GB
# CPUs: 8

set -eo pipefail

echo "=========================================="
echo "Workflow 3: Jwalk Collection"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo ""

# Navigate to project directory
cd /storage/group/epo2/default/ims86/git_repos/EntDetect

# Activate environment
source ~/.bashrc
conda activate entdetect

# Configuration
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3
XP_DIR=$DATASTORE/outputs/workflow2/OP_AA/XP

# Verify prerequisites
echo "Verifying prerequisites..."
if [[ ! -d "$XP_DIR" ]]; then
    echo "ERROR: XP directory not found at $XP_DIR"
    exit 1
fi

n_xp_files=$(ls -1 $XP_DIR/*_XP 2>/dev/null | wc -l)
if [[ $n_xp_files -eq 0 ]]; then
    echo "ERROR: No XP files found in $XP_DIR"
    exit 1
fi
echo "✓ Found $n_xp_files XP files"
echo ""

# Create output directory
mkdir -p $OUTDIR/logs

echo "Collecting Jwalk/XP outputs..."
echo "XP directory: $XP_DIR"
echo "Output directory: $OUTDIR"
echo ""

# Run collection
python scripts/collect_jwalk.py \
    --xp_dir "$XP_DIR" \
    --n_traj 1000 \
    --n_frames 335 \
    --outdir "$OUTDIR" \
    --ID 1ZMR \
    --prot_len 387 \
    2>&1 | tee "$OUTDIR/logs/jwalk_collection.log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=========================================="
echo "Completion time: $(date)"
echo "=========================================="

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✓ Jwalk collection completed successfully!"
    echo ""
    echo "Next step: Submit consistency test job"
    echo "  sbatch assets/slurm/scripts/run_workflow3_consistency.slurm"
else
    echo "✗ Jwalk collection failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE

