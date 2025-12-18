from EntDetect.statistics import ProteomeLogisticRegression

"""
python scripts/run_population_modeling.py 
--dataframe_files /path/to/dataframe/files/ 
--outdir TestingGrounds/ProteomeLogisticRegression/ 
--gene_list /path/to/gene_list.txt 
--tag test_regression 
--reg_formula "cut_C_Rall ~ AA + region"
"""

if __name__ == "__main__":

    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--dataframe_files", type=str, required=True, help="Path to the directory containing the dataframe files")
    parser.add_argument("--outdir", type=str, required=True, help="Path to the output directory for regression results")
    parser.add_argument("--gene_list", type=str, required=True, help="Path to the gene list file")
    parser.add_argument("--tag", type=str, required=True, help="Tag for the regression results")
    parser.add_argument("--reg_formula", type=str, default='cut_C_Rall ~ AA + region', help="Regression formula to use")
    args = parser.parse_args()
    print(args)
    dataframe_files = args.dataframe_files
    outdir = args.outdir
    gene_list = args.gene_list
    tag = args.tag
    reg_formula = args.reg_formula

    ## initialize the regression object
    ProtRegession = ProteomeLogisticRegression(dataframe_files=dataframe_files,
                                                outdir=outdir,
                                                gene_list=gene_list,
                                                tag=tag,
                                                reg_formula=reg_formula)
    print(ProtRegession)

    ## do the clustering
    ProtRegession.load_data(sep='|', reg_var=['AA', 'region'], response_var='cut_C_Rall', var2binarize=['cut_C_Rall', 'region'], mask_column='mapped_resid')

    ## Run the regression analysis
    reg_df = ProtRegession.run()

    # Define output files and get gene list
    reg_outfile = os.path.join(outdir, f"regression_results_{tag}.csv")
    reg_df.to_csv(reg_outfile, index=False, sep='|')
    print(f"Regression results saved to {reg_outfile}")
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
