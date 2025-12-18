from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements, MSMNonNativeEntanglementClustering
from EntDetect.order_params import CalculateOP

"""
python scripts/run_MSM.py 
--outdir TestingGrounds/MSM/ 
--ID 1ZMR_prod 
--OPpath /path/to/OP/data/ 
--start 0 
--n_large_states 10 
--lagtime 20 
--rm_traj_list 65 75 155 162
"""

if __name__ == "__main__":

    import multiprocessing as mp
    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outdir", type=str, required=True, help="Path to output directory")
    parser.add_argument("--OPpath", type=str, required=True, help="Path to directory containing G and Q directories created by GQ.py")
    parser.add_argument("--ID", type=str, required=True, help="base name for output files")
    parser.add_argument("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
    parser.add_argument("--n_large_states", type=int, required=False, help="Number of large states for MSM", default=10)
    parser.add_argument("--lagtime", type=int, required=False, help="Lag time for MSM", default=20)
    parser.add_argument("--rm_traj_list", type=int, nargs='+', required=False, help="List of trajectory indices to remove", default=[])
    args = parser.parse_args()
    print(args)
    outdir = args.outdir
    OPpath = args.OPpath
    ID = args.ID
    start = args.start
    n_large_states = args.n_large_states
    lagtime = args.lagtime
    rm_traj_list = args.rm_traj_list

    ## Load the MSM object with the parameters
    MSM = MSMNonNativeEntanglementClustering(outdir=outdir,
                         ID=ID,
                         OPpath=OPpath,
                         start=start, 
                         n_large_states=n_large_states,
                         rm_traj_list=rm_traj_list,
                         lagtime=lagtime)
    print(MSM)
    MSM.run()

    ## Sample representative structures from the MSM
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
