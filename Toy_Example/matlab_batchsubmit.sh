#!/bin/bash
#SBATCH -J sOFMM_latent_scen6_14	# Job name for the array
#SBATCH -o sOFMM_latent_scen6_14_%A.out  # Shared standard output with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 1300	# Memory request
#SBATCH -t 0-05:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
scenarios=(6 14)
iter=1
for scenario in ${scenarios[@]}; do
	echo This is scenario ${scenario} iteration ${iter} sample ${SLURM_ARRAY_TASK_ID}.$'\n'
	matlab -nodisplay -nosplash -r "sOFMM_main_latent(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
done
