#! /bin/bash
#SBATCH -J sim_SRS 	      # Job name
#SBATCH -o sim_SRS_%A.out  # Shared std out for job array
#SBATCH -p shared             # Partition
#SBATCH -n 1                  # Number of nodes
#SBATCH --mem 1000             # Memory request in MB
#SBATCH -t 0-00:05:00         # Runtime (D-HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Email account

module load matlab
samp_n=1
for scenario in $(seq 1 4); do
	echo This is population scenario ${scenario} iteration ${SLURM_ARRAY_TASK_ID} sample ${samp_n}.$'\n'
matlab -nodisplay -nosplash -r "sample_SRS(${scenario},${SLURM_ARRAY_TASK_ID},${samp_n})" 
