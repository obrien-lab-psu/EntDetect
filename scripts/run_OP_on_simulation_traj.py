from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements
from EntDetect.order_params import CalculateOP

"""
python scripts/run_OP_on_simulation_traj.py 
--Traj 1 
--PSF /storage/group/epo2/default/for_Ian_from_Yang/1zmr_model_clean_ca.psf 
--DCD /storage/group/epo2/default/for_Ian_from_Yang/prod_traj_reduce_saving/1_prod.dcd 
--ID 1ZMR 
--COR /storage/group/epo2/default/for_Ian_from_Yang/1zmr_model_clean_ca.cor 
--sec_elements /storage/group/epo2/default/for_Ian_from_Yang/secondary_struc_defs.txt 
--domain /storage/group/epo2/default/for_Ian_from_Yang/domain_def.dat 
--outdir TestingGrounds/run_OP_on_simulation_traj/ 
--start 6660
"""

if __name__ == "__main__":

    import multiprocessing as mp
    import sys,os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--Traj", type=str, required=True, help="Trajectory number to use")
    parser.add_argument("--PSF", type=str, required=True, help="Path to to reference structure PSF file")
    parser.add_argument("--DCD", type=str, required=True, help="Path to DCD file of trajectory")
    parser.add_argument("--ID", type=str, required=True, help="An id for the PDB to be analyzed")
    parser.add_argument("--COR", type=str, required=True, help="Path to reference structure COR file")
    parser.add_argument("--sec_elements", type=str, required=True, help="Path to secondary structure definitions file")
    parser.add_argument("--domain", type=str, required=True, help="Path to domain definitions file")
    parser.add_argument("--start", type=int, required=False, help="Frame number to start analysis from", default=0)
    parser.add_argument("--outdir", type=str, required=False, help="Output directory for results", default='./')    
    args = parser.parse_args()
    print(args)
    traj = args.Traj
    ID = args.ID
    psf = args.PSF
    dcd = args.DCD
    cor = args.COR
    sec_elements = args.sec_elements
    domain = args.domain
    start = args.start
    outdir = args.outdir

    ## Load the Order Parameter object with the file parameters
    CalcOP = CalculateOP(outdir=outdir,
                         Traj=traj,
                         ID=ID,
                         psf=psf,
                         cor=cor,
                         sec_elements=sec_elements,
                         dcd=dcd,
                         domain=domain,
                         start=start) #
    print(CalcOP)

    ## Calculate the fraction of native contacts (Q)
    Qdata_dict = CalcOP.Q()
    print(Qdata_dict.keys())

    ## Calculate the fraction of native contacts with a change of entanglement (G) and all associated entanglement features 
    # Case #1: 
    #   - Use topoly to find crossings and use the TLN to assess entanglement changes
    #   - Use the alpha carbons to define native contacts
    #   - Trajectory is coarse-grained (CG) 
    Gdata_dict = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10) 
    print(Gdata_dict.keys())

    # Case #2: 
    #   - Skip topoly to find crossings and use the GLN to assess entanglement changes
    #   - Use the alpha carbons to define native contacts
    #   - Trajectory is coarse-grained (CG) 
    # Gdata_dict = CalcOP.G(topoly=False, Calpha=True, CG=True, nproc=10) 
    # print(Gdata_dict.keys())

    # Case #3: 
    #   - Skip topoly to find crossings and use the GLN to assess entanglement changes
    #   - Use all heavy atoms to define native contacts
    #   - Trajectory is all-atom (not CG) 
    # Gdata_dict = CalcOP.G(topoly=False, Calpha=False, CG=False, nproc=10) 
    # print(Gdata_dict.keys())

    # Case #4: 
    #   - Use topoly to find crossings and use the TLN to assess entanglement changes
    #   - Use all heavy atoms to define native contacts
    #   - Trajectory is all-atom (not CG) 
    # Gdata_dict = CalcOP.G(topoly=True, Calpha=False, CG=False, nproc=10) 
    # print(Gdata_dict.keys())

    ## Calculate the mirror symmetry order parameter K
    Kdata_dict = CalcOP.K()
    print(Kdata_dict.keys())

    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
