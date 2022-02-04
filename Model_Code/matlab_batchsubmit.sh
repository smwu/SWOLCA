#!/bin/bash
#SBATCH -J wsOFMM_iter1	# Job name for the array
#SBATCH -o wsOFMM_iter1_%A.out  # Shared standard output with job ID
#SBATCH -e wsOFMM_iter1_%A.err  # Shared standard error with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 1300	# Memory request
#SBATCH -t 0-7:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
iter=1
for scenario in $(seq 5 16); do
	echo This is scenario ${scenario} iteration ${iter} sample ${SLURM_ARRAY_TASK_ID}.$'\n'
	matlab -nodisplay -nosplash -r "wsOFMM_main(${scenario},${iter},${SLURM_ARRAY_TASK_ID})"
done
