#!/bin/bash
#SBATCH -J Uneq_scen201_202_wOFMM	# Job name for the array
#SBATCH -o Uneq_scen201_202_wOFMM%A.out  # Shared standard output with job ID
#SBATCH -e Uneq_scen201_202_wOFMM%A.err  # Shared standard error with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 2500	# Memory request
#SBATCH -t 0-04:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab/R2022a-fasrc01
scenarios=(201 202)
iter=1
for scenario in ${scenarios[@]}; do
	echo This is scenario ${scenario} iteration ${iter} sample ${SLURM_ARRAY_TASK_ID}.$'\n'
	# matlab -nodisplay -nosplash -r "sOFMM_main(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
	# matlab -nodisplay -nosplash -r "sOFMM_main_latent(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
        # matlab -nodisplay -nosplash -r "wsOFMM_main(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
        matlab -nodisplay -nosplash -r "wOFMM_main_latent(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
done
