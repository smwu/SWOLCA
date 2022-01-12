#!/bin/bash
#SBATCH -J wsOFMM	# Job name for the array
#SBATCH -o wsOFMM_%a.out  # Standard output with job array index
#SBATCH -e wsOFMM_%a.err  # Standard error with job array index
#SBATCH -p shared       # Partition to submit to
#SBATCH -n 1		# Number of cores
#SBATCH --mem 12000	# Memory request (12Gb)
#SBATCH -t 2-8:00:00	# Runtime (D-HH:MM) days, hours, mins
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications 
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load matlab
matlab  -nodisplay -nosplash < wsOFMM"${SLURM_ARRAY_TASK_ID}".m
