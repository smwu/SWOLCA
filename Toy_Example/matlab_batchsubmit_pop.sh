#!/bin/bash
#SBATCH -J sOFMM_latent_scen2	# Job name for the array
#SBATCH -o sOFMM_latent_scen2_%A.out  # Shared standard output with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 8000	# Memory request
#SBATCH -t 0-05:0:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
samp_n=0
for scenario in $(seq 2 2); do
	echo This is scenario ${scenario} iteration ${SLURM_ARRAY_TASK_ID}.$'\n'
	matlab -nodisplay -nosplash -r "sOFMM_main_latent(${scenario},${SLURM_ARRAY_TASK_ID},${samp_n})"
done
