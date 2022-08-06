#! /bin/bash
#SBATCH -J sim_uneq_pi_subpop 	      # Job name
#SBATCH -o sim_uneq_pi_subpop_%A.out  # Shared std out for job array
#SBATCH -p shared             # Partition
#SBATCH -n 1                  # Number of nodes
#SBATCH --mem 2000             # Memory request in MB
#SBATCH -t 0-00:01:00         # Runtime (D-HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Email account

module load matlab
echo This is dataset ${SLURM_ARRAY_TASK_ID}.$'\n'
matlab -nodisplay -nosplash -r "sim_uneq_pi_subpop(${SLURM_ARRAY_TASK_ID})" 
