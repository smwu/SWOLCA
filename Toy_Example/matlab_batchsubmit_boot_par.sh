#!/bin/bash
#SBATCH -o B_M200_SRSWOR_2%A.out
#SBATCH -e B_M200_SRSWOR_2%A.err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem 18000
#SBATCH -t 1-10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=stephaniewu@fas.harvard.edu

module load matlab
for rep in {51..100}; do
	matlab  -nodisplay -nosplash -r "bootPRS_sim_par($rep,${SLURM_ARRAY_TASK_ID}, 'uOFMsimdata_B', 'bPRS_simResultsB')"
done
