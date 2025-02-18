#!/bin/bash
#SBATCH --job-name=openmm_sim
#SBATCH --output=simulation_%j.out
#SBATCH --error=simulation_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# Load required modules (modify according to your HPC setup)
module load cuda/11.4
module load anaconda3

# Activate conda environment
source activate openmm-dev

# Run simulation
python scripts/local_run.py --pdb $1 --steps $2
