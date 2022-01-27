#! /bin/bash
#SBATCH -J sim_equal 	      # Job name
#SBATCH -o sim_equal%A.out    # Shared std out for job array
#SBATCH -p shared             # Partition
#SBATCH -n 1                  # Number of nodes
#SBATCH --mem 800             # Memory request in MB
#SBATCH -t 0-00:05:00         # Runtime (D-HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Email account

module load matlab
echo This is dataset ${SLURM_ARRAY_TASK_ID}.
matlab -nodisplay -nosplash -r "sim_equal_subpop(${SLURM_ARRAY_TASK_ID})" 
