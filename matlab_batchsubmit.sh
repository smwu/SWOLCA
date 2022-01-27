#!/bin/bash
#SBATCH -J wsOFMM	# Job name for the array
#SBATCH -o wsOFMM_%A.out  # Shared standard output with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 2000	# Memory request
#SBATCH -t 0-15:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
for scenario in $(seq 2 4); do
	echo This is scenario ${scenario} iteration ${SLURM_ARRAY_TASK_ID}.$'\n'
	matlab -nodisplay -nosplash -r "wsOFMM_main(${scenario},${SLURM_ARRAY_TASK_ID})"
done
