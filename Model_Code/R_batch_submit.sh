#!/bin/bash
#SBATCH -J Rcpp_WOLCA_scen111211_        # Job name for the array
#SBATCH -o Rcpp_WOLCA_scen111211_%A.out  # Shared standard output with job ID
#SBATCH -e Rcpp_WOLCA_scen111211_%A.err  # Shared standard error with job ID
#SBATCH -p shared      # Partition to submit to
#SBATCH -n 1	       # Number of cores
#SBATCH -t 0-8:00:00  # Runtime (D-HH:MM:SS)
#SBATCH --mem=5000     # Memory request
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load R/4.2.2-fasrc01 gcc/10.2.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
scen_samp=111211
iter_pop=1
covs=2
Rscript /n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/Model_Code/WOLCA_main_Rcpp.R ${scen_samp} ${iter_pop} ${SLURM_ARRAY_TASK_ID} ${covs}
