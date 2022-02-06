#!/bin/bash
#SBATCH -J wsOFMM_scen5_16	# Job name for the array
#SBATCH -o wsOFMM_scen5_16_%A.out  # Shared standard output with job ID
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 1300	# Memory request
#SBATCH -t 0-8:00:00	# Runtime (D-HH:MM:SS) days, hours, mins, secs
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
samp_n=1
for scenario in $(seq 5 16); do
	echo This is scenario ${scenario} iteration ${SLURM_ARRAY_TASK_ID} sample ${samp_n}.$'\n'
	matlab -nodisplay -nosplash -r "wsOFMM_main(${scenario},${SLURM_ARRAY_TASK_ID},${samp_n})"
done
