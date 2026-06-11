#!/usr/bin/env python
"""
Collect per-trajectory Jwalk (XP) outputs into a consolidated Jwalk.npy array.

This is typically run as a separate SLURM job (takes ~3-4 hours for 1000 trajectories).

Usage:
    python scripts/collect_jwalk.py \
        --xp_dir /path/to/OP_AA/XP \
        --n_traj 1000 \
        --n_frames 335 \
        --outdir /path/to/output/ \
        --ID 1ZMR \
        --prot_len 387
"""

import argparse
import logging
from pathlib import Path
from EntDetect.order_params import CollectOP

def main(argv=None):
    parser = argparse.ArgumentParser(description="Collect per-trajectory Jwalk/XP files")
    parser.add_argument("--xp_dir", type=str, required=True, help="Directory of per-traj XP files")
    parser.add_argument("--n_traj", type=int, required=True, help="Total number of trajectories")
    parser.add_argument("--n_frames", type=int, required=True, help="Frames per trajectory")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--ID", type=str, default="1ZMR", help="ID prefix for trajectory files")
    parser.add_argument("--prot_len", type=int, required=True, help="Protein length")
    args = parser.parse_args(argv)
    
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s'
    )
    log = logging.getLogger(__name__)
    
    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    log.info(f"Collecting Jwalk from {args.xp_dir}")
    log.info(f"  Trajectories: {args.n_traj}")
    log.info(f"  Frames per traj: {args.n_frames}")
    log.info(f"  Protein length: {args.prot_len}")
    log.info(f"  Output: {outdir}")
    
    # Initialize collector
    collector = CollectOP(
        sasa_dir=None,
        xp_dir=args.xp_dir,
        outdir=str(outdir),
        ID=args.ID,
        n_traj=args.n_traj,
        n_frames=args.n_frames,
        prot_len=args.prot_len,
    )
    
    # Collect Jwalk
    log.info("Starting Jwalk collection...")
    jwalk_file = collector.collect_Jwalk()
    log.info(f"✓ Jwalk collection complete: {jwalk_file}")
    
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
