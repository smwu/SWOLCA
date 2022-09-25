#!/bin/bash
#SBATCH -J Uneq_scen201_202_bwsOFMM_P1
#SBATCH -o Uneq_scen201_202_bwsOFMM_P1%A.out
#SBATCH -e Uneq_scen201_202_bwsOFMM_P1%A.err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem 18000
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=stephaniewu@fas.harvard.edu

module load matlab
scenarios=(201 202)
iter=1
for scenario in ${scenarios[@]}; do
	for samp_n in {1..25}; do
		matlab  -nodisplay -nosplash -r "bwsOFMM_main_latent(${scenario},${iter},${samp_n},${SLURM_ARRAY_TASK_ID})"
	done
done
