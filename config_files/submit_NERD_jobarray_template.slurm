#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=nerdYY
#SBATCH --partition inf
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G    # Required RAM per cpu (default is 4 Gb)
#SBATCH --time=3-1:30			  # Max runtime (d-h:mm) of task
#SBATCH --output=nerd_XXXm_YYpx_%A-%a.out	  # %A jobID, %a array task
#SBATCH --error=nerd_XXXm_YYpx_%A-%a.err
#SBATCH --array=1-137%15       # set array numbers, but limit the number of simultaneious running tasks to 15

# set python environment
module load Python
source /net/labdata/maaike/Environments/NERD/bin/activate

# find all image files that should be processed
FILES=(/net/labdata/maaike/Data/S1_SAR/relorbs/2021-22/*.tif)

# output jobArr ID to .out for administration
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# get config file
config_file=/net/labdata/maaike/NERD/python/config_files/config_NERD_image_template.ini

# apply NERD to all images
srun python -u /net/labdata/maaike/NERD/python/scripts/main.py $config_file ${FILES[$SLURM_ARRAY_TASK_ID]}
