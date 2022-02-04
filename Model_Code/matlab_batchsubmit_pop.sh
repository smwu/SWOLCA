#!/bin/bash
#SBATCH -J wsOFMM_pop	# Job name for the array
#SBATCH -o wsOFMM_pop%A.out  # Shared standard output with job ID
#SBATCH -e wsOFMM_pop%A.err  # Shared standard error with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 10000	# Memory request
#SBATCH -t 0-15:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
samp_n=0
for scenario in $(seq 1 4); do
	echo This is scenario ${scenario} iteration ${SLURM_ARRAY_TASK_ID}.$'\n'
	matlab -nodisplay -nosplash -r "wsOFMM_main(${scenario},${SLURM_ARRAY_TASK_ID},${samp_n})"
done
