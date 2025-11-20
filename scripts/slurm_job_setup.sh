#!/bin/bash -l
# Job Name:
#SBATCH -J python
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3988
# Wall clock limit:
#SBATCH --time=02:00:00

module purge
module load GCC/11.3.0 Python/3.10.4 IPython/8.5.0 OpenMPI/4.1.4 astropy/5.1.1 h5py/3.7.0 matplotlib/3.5.2 tqdm

# Run the python script and give the script the index of this specific job
python setup_slurm.py ${SLURM_ARRAY_TASK_ID} ${SLURM_CPUS_PER_TASK}

