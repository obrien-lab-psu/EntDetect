from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements, ClusterNonNativeEntanglements

"""
python scripts/run_nonnative_entanglement_clustering.py 
--outdir TestingGrounds/ClusterNonNativeEntanglements/ 
--pkl_file_path /path/to/pkl/files/ 
--trajnum2pklfile_path /path/to/mapping.txt 
--traj_dir_prefix /path/to/traj/prefix
"""

if __name__ == "__main__":

    import multiprocessing as mp
    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outdir", type=str, required=False, help="Output directory for results", default='./')
    parser.add_argument("--pkl_file_path", type=str, required=False, help="Path to directory containing pickled non-native entanglement data", default=None)
    parser.add_argument("--trajnum2pklfile_path", type=str, required=False, help="Path to file mapping trajectory numbers to pkl files", default=None)
    parser.add_argument("--traj_dir_prefix", type=str, required=False, help="Path prefix to directory containing trajectories", default=None)
    args = parser.parse_args()
    print(args)
    outdir = args.outdir
    pkl_file_path = args.pkl_file_path
    trajnum2pklfile_path = args.trajnum2pklfile_path
    traj_dir_prefix = args.traj_dir_prefix

    ## initialize the clustering object
    clustering_NNents = ClusterNonNativeEntanglements(pkl_file_path=pkl_file_path, 
                                                      trajnum2pklfile_path=trajnum2pklfile_path,
                                                      traj_dir_prefix=traj_dir_prefix,
                                                      outdir=outdir)
  
    ## do the clustering
    clustering_NNents.cluster(start_frame=6600)
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
