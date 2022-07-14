#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=expGEE		# the name of your Job as it appears in the queue.
#SBATCH --partition inf
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G			# Required RAM for your job.
#SBATCH --time=3-1:30			# Max runtime (d-h:mm) of task
#SBATCH --output=output_GEE_export%A.out	#  %A will be replaced by your jobID.
#SBATCH --error=error_GEE_export%A.err		#  %A will be replaced by your jobID

module load Python
source /net/labdata/maaike/Environments/NERD/bin/activate

python -u /net/labdata/maaike/NERD/python/scripts/gee_to_python.py
