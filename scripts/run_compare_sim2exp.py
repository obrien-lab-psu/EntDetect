from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements, MSMNonNativeEntanglementClustering
from EntDetect.order_params import CalculateOP
from EntDetect.compare_sim2exp import MassSpec

"""
python scripts/run_compare_sim2exp.py 
--msm_data_file /path/to/msm_data.csv 
--meta_dist_file /path/to/meta_dist.csv 
--LiPMS_exp_file /path/to/LiPMS_exp.csv 
--sasa_data_file /path/to/sasa_data.csv 
--XLMS_exp_file /path/to/XLMS_exp.csv 
--dist_data_file /path/to/dist_data.csv 
--cluster_data_file /path/to/cluster_data.csv 
--OPpath /path/to/OP/ 
--AAdcd_dir /path/to/AAdcd/ 
--native_AA_pdb /path/to/native.pdb 
--state_idx_list 1 2 3 4 5 
--prot_len 390 
--last_num_frames 100 
--rm_traj_list 65 75 155 
--native_state_idx 0 
--outdir TestingGrounds/Compare/ 
--ID 1ZMR 
--start 0 
--end 1000 
--stride 1 
--num_perm 1000 
--n_boot 100 
--lag_frame 20 
--nproc 10
"""

if __name__ == "__main__":

    import multiprocessing as mp
    import sys, os
    import argparse
    import time

    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--msm_data_file", type=str, required=True, help="Path to MSM mapping file")
    parser.add_argument("--meta_dist_file", type=str, required=True, help="Path to meta-distance file")
    parser.add_argument("--LiPMS_exp_file", type=str, required=True, help="Path to LiP-MS experimental data file")
    parser.add_argument("--sasa_data_file", type=str, required=True, help="Path to SASA data file")
    parser.add_argument("--XLMS_exp_file", type=str, required=True, help="Path to XL-MS experimental data file")
    parser.add_argument("--dist_data_file", type=str, required=True, help="Path to distance data file")
    parser.add_argument("--cluster_data_file", type=str, required=True, help="Path to clustering data file")
    parser.add_argument("--OPpath", type=str, required=True, help="Path to order parameters directory")
    parser.add_argument("--AAdcd_dir", type=str, required=True, help="Path to all-atom DCD files directory")
    parser.add_argument("--native_AA_pdb", type=str, required=True, help="Path to native all-atom PDB file")
    parser.add_argument("--state_idx_list", type=int, nargs='+', required=True, help="List of state indices to analyze")
    parser.add_argument("--prot_len", type=int, required=True, help="Length of the protein")
    parser.add_argument("--last_num_frames", type=int, required=True, help="Number of last frames to consider")
    parser.add_argument("--rm_traj_list", type=int, nargs='+', required=True, help="List of trajectory indices to remove")
    parser.add_argument("--native_state_idx", type=int, required=True, help="Index of the native state")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory for results")
    parser.add_argument("--ID", type=str, required=True, help="An ID for the analysis")
    parser.add_argument("--start", type=int, required=True, help="Start frame index")
    parser.add_argument("--end", type=int, required=True, help="End frame index")
    parser.add_argument("--stride", type=int, required=True, help="Stride for frame selection")
    parser.add_argument("--verbose", action='store_true', help="Enable verbose output")
    parser.add_argument("--num_perm", type=int, required=True, help="Number of permutations for statistical tests")
    parser.add_argument("--n_boot", type=int, required=True, help="Number of bootstrap samples")
    parser.add_argument("--lag_frame", type=int, required=True, help="Lag time in frames")
    parser.add_argument("--nproc", type=int, required=True, help="Number of processes for parallel computation")
    args = parser.parse_args()
    print(args)
    msm_data_file = args.msm_data_file
    meta_dist_file = args.meta_dist_file
    LiPMS_exp_file = args.LiPMS_exp_file
    sasa_data_file = args.sasa_data_file
    XLMS_exp_file = args.XLMS_exp_file
    dist_data_file = args.dist_data_file
    cluster_data_file = args.cluster_data_file
    OPpath = args.OPpath
    AAdcd_dir = args.AAdcd_dir
    native_AA_pdb = args.native_AA_pdb
    state_idx_list = args.state_idx_list
    prot_len = args.prot_len
    last_num_frames = args.last_num_frames
    rm_traj_list = args.rm_traj_list
    native_state_idx = args.native_state_idx
    outdir = args.outdir
    ID = args.ID
    start = args.start
    end = args.end
    stride = args.stride
    verbose = args.verbose
    num_perm = args.num_perm
    n_boot = args.n_boot
    lag_frame = args.lag_frame
    nproc = args.nproc

    ## initialize the MassSpec object
    MS = MassSpec(msm_data_file=msm_data_file,
                    meta_dist_file=meta_dist_file,
                    LiPMS_exp_file=LiPMS_exp_file,
                    sasa_data_file=sasa_data_file,
                    XLMS_exp_file=XLMS_exp_file,
                    dist_data_file=dist_data_file,
                    cluster_data_file=cluster_data_file,
                    OPpath=OPpath,
                    AAdcd_dir=AAdcd_dir,
                    native_AA_pdb=native_AA_pdb,
                    state_idx_list=state_idx_list,
                    prot_len=prot_len,
                    last_num_frames=last_num_frames,
                    rm_traj_list=rm_traj_list,
                    native_state_idx=native_state_idx,
                    outdir=outdir,
                    ID=ID,
                    start=start,
                    end=end,
                    stride=stride,
                    verbose=verbose,
                    num_perm=num_perm,
                    n_boot=n_boot,
                    lag_frame=lag_frame,
                    nproc=nproc)

    # run the consistency test
    consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
    print(f'consist_data_file: {consist_data_file}')
    print(f'consist_result_file: {consist_result_file}')

    # select the representative structures from the consistency test
    MS.select_rep_structs(consist_data_file, consist_result_file, total_traj_num_frames=335, last_num_frames=67)
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
