from EntDetect.statistics import MonteCarlo

"""
python scripts/run_montecarlo.py 
--dataframe_files /path/to/dataframe/files/ 
--outpath TestingGrounds/MonteCarlo/ 
--gene_list /path/to/gene_list.txt 
--tag test_run 
--steps 100000 
--n_groups 4 
--C1 1.0 
--C2 2.5 
--beta 0.05
"""

if __name__ == "__main__":

    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--dataframe_files", type=str, required=True, help="Path to the directory containing the dataframe files")
    parser.add_argument("--outpath", type=str, required=True, help="Path to the output directory for Monte Carlo results")
    parser.add_argument("--gene_list", type=str, required=True, help="Path to the gene list file")
    parser.add_argument("--tag", type=str, required=True, help="Tag for the Monte Carlo results")
    parser.add_argument("--reg_formula", type=str, default='cut_C_Rall ~ region + AA', help="Regression formula to use")
    parser.add_argument("--response_var", type=str, default='cut_C_Rall', help="Response variable for the regression")
    parser.add_argument("--test_var", type=str, default='region', help="Variable to test in the regression")
    parser.add_argument("--random", action='store_true', help="Use random sampling for the Monte Carlo simulation")
    parser.add_argument("--n_groups", type=int, default=4, help="Number of groups for the Monte Carlo simulation")
    parser.add_argument("--steps", type=int, default=100000, help="Number of steps for the Monte Carlo simulation")
    parser.add_argument("--C1", type=float, default=1.0, help="C1 parameter for the Monte Carlo simulation")
    parser.add_argument("--C2", type=float, default=2.5, help="C2 parameter for the Monte Carlo simulation")
    parser.add_argument("--beta", type=float, default=0.05, help="Beta parameter for the Monte Carlo simulation")
    parser.add_argument("--linearT", action='store_true', help="Use linear temperature for the Monte Carlo simulation")
    args = parser.parse_args()
    print(args)
    dataframe_files = args.dataframe_files
    outpath = args.outpath
    gene_list = args.gene_list
    tag = args.tag
    reg_formula = args.reg_formula
    response_var = args.response_var
    test_var = args.test_var
    random = args.random
    n_groups = args.n_groups
    steps = args.steps
    C1 = args.C1
    C2 = args.C2
    beta = args.beta
    linearT = args.linearT

    ## initialize the clustering object
    MC = MonteCarlo(dataframe_files=dataframe_files, 
                    outpath=outpath, 
                    gene_list=gene_list, 
                    tag=tag, 
                    reg_formula=reg_formula,
                    response_var=response_var, 
                    test_var=test_var,
                    random=random, 
                    n_groups=n_groups, 
                    steps=steps, 
                    C1=C1, 
                    C2=C2, 
                    beta=beta,
                    linearT=linearT)
    print(MC)

    ## Load the data into the MonteCarlo object
    MC.load_data(sep='|', reg_var=['AA', 'region'], response_var='cut_C_Rall', var2binarize=['cut_C_Rall', 'region'], mask_column='mapped_resid', ID_column='gene', Length_column='uniprot_length')

    ## Run the Monte Carlo simulation
    MC.run(encoded_df=MC.data, ID_column='gene')
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
