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

def main(argv=None):

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
    parser.add_argument("--ent_detection_method", type=int, required=False, help="ENT detection method: 1=any GLN, 2=any TLN (default), 3=both GLN and TLN same termini", default=2)
    parser.add_argument(
        "--resolution",
        type=str,
        choices=["cg", "aa"],
        default="cg",
        help="Trajectory resolution: 'cg' (coarse-grained, default) or 'aa' (all-atom)",
    )
    parser.add_argument(
        "--contacts",
        type=str,
        choices=["calpha", "heavy"],
        default=None,
        help="Contact definition to use for OP/G: 'calpha' or 'heavy'. Default: calpha for cg, heavy for aa.",
    )
    parser.add_argument(
        "--no_topoly",
        action="store_true",
        help="Disable topoly crossing detection (uses GLN-only workflow). Default: topoly enabled.",
    )
    parser.add_argument("--nproc", type=int, required=False, default=10, help="Number of processes for G calculation (default: 10)")
    parser.add_argument("--outdir", type=str, required=False, help="Output directory for results", default='./')    
    args = parser.parse_args(argv)
    print(args)
    traj = args.Traj
    ID = args.ID
    psf = args.PSF
    dcd = args.DCD
    cor = args.COR
    sec_elements = args.sec_elements
    domain = args.domain
    start = args.start
    ent_detection_method = args.ent_detection_method
    outdir = args.outdir

    if args.contacts is None:
        contacts = "calpha" if args.resolution == "cg" else "heavy"
    else:
        contacts = args.contacts

    Calpha = contacts == "calpha"
    CG = args.resolution == "cg"
    topoly = not args.no_topoly
    nproc = args.nproc

    ## Load the Order Parameter object with the file parameters
    CalcOP = CalculateOP(outdir=outdir,
                         Traj=traj,
                         ID=ID,
                         psf=psf,
                         cor=cor,
                         sec_elements=sec_elements,
                         dcd=dcd,
                         domain=domain,
                         start=start,
                         ent_detection_method=ent_detection_method) #
    print(CalcOP)

    ## Calculate the fraction of native contacts (Q)
    Qdata_dict = CalcOP.Q()
    print(Qdata_dict.keys())

    ## Calculate the fraction of native contacts with a change of entanglement (G) and all associated entanglement features
    # Defaults preserve historical behavior: coarse-grained trajectory, Calpha contacts, topoly enabled.
    Gdata_dict = CalcOP.G(topoly=topoly, Calpha=Calpha, CG=CG, nproc=nproc)
    print(Gdata_dict.keys())

    ## Calculate the mirror symmetry order parameter K
    Kdata_dict = CalcOP.K()
    print(Kdata_dict.keys())

    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
