from EntDetect.statistics import FoldingPathwayStats
import pandas as pd

"""
python scripts/run_Foldingpathway.py 
--msm_data_file MSMfull/1ZMR_prod_meta_set_A80%Native.csv 
--meta_set_file MSMfull/1ZMR_prod_meta_set.csv 
--traj_type_col traj_type_A80%Native 
--outdir Foldingpathway_A80%Native 
--rm_traj_list 65 75 155 162
"""

if __name__ == "__main__":

    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--msm_data_file", type=str, required=True, help="Path to MSM data file")
    parser.add_argument("--meta_set_file", type=str, required=True, help="Path to meta set file")
    parser.add_argument("--traj_type_col", type=str, required=True, help="Trajectory type column name")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory for results")
    parser.add_argument("--rm_traj_list", type=int, nargs='+', required=False, help="List of trajectory indices to remove", default=[])
    args = parser.parse_args()
    print(args)
    msm_data_file = args.msm_data_file
    meta_set_file = args.meta_set_file
    traj_type_col = args.traj_type_col
    outdir = args.outdir
    rm_traj_list = args.rm_traj_list

    # Load MSM data
    msm_data = pd.read_csv(msm_data_file)

    ## initialize the folding pathway stats object
    FP = FoldingPathwayStats(msm_data=msm_data,
                             meta_set_file=meta_set_file,
                             tarj_type_col=traj_type_col,
                             outdir=outdir,
                             rm_traj_list=rm_traj_list)

    ## get the post-transitional folding pathways
    folding_pathways = FP.post_trans()

    ## JS divergence
    JS_divergence = FP.JS_divergence()

    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
