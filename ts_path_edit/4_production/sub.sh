#!/bin/bash
#SBATCH -A <account>
#SBATCH --output=time.out
#SBATCH --error=error.out
#SBATCH --job-name=p_mech_gen
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --time=04:00:00

## Perlmutter script example
cd $SLURM_SUBMIT_DIR

module load python
conda activate <path_to_conda>/conda_env

python swap_production.py --nprocs 128
