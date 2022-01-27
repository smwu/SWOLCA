#! /bin/bash
#SBATCH -J parsim_pi_subpop   # Job name
#SBATCH -o parsim_pi.out      # Std out
#SBATCH -e parsim_pi.err      # Std err
#SBATCH -p shared                # Partition
#SBATCH -c 32		      	 # Number of cores
#SBATCH -n 1                  	 # Number of nodes
#SBATCH --mem 12000           	 # Memory request (12Gb)
#SBATCH -t 7-00:00            	 # Runtime (D-HH:MM)
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Email account

module load matlab
matlab -nodisplay -nosplash < parsim_wsRPC_pi_subpop.m 
