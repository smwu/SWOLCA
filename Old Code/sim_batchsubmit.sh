#! /bin/bash
#SBATCH -J sim_equal_subpop   # Job name
#SBATCH -o sim_equal.out      # Std out
#SBATCH -e sim_equal.err      # Std err
#SBATCH -p shared             # Partition
#SBATCH -n 1                  # Number of nodes
#SBATCH --mem 12000           # Memory request (12Gb)
#SBATCH -t 2-00:00            # Runtime (D-HH:MM)
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Email account

module load matlab
matlab -nodisplay -nosplash < sim_wsRPC_equal_subpop.m 
