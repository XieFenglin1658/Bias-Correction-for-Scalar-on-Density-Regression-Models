#!/bin/bash
#SBATCH --account=5sigma
#SBATCH --job-name=nhanes_simex_scaled
#SBATCH --output=logs/simex_%A_%a.out
#SBATCH --error=logs/simex_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-70

cd /insomnia001/depts/5sigma/users/fx2212/simulation_nhanes
mkdir -p logs

source /insomnia001/home/fx2212/miniconda3/bin/activate fda-43

M0_VALUES=(10 20 50 100 200 500 1000)
N_BLOCKS=10
BLOCK_SIZE=10

idx=$((SLURM_ARRAY_TASK_ID - 1))
m0_idx=$((idx / N_BLOCKS))
block_idx=$((idx % N_BLOCKS))

m0=${M0_VALUES[$m0_idx]}
start_rep=$((block_idx * BLOCK_SIZE + 1))
B_REPS=$BLOCK_SIZE

echo "Processing m0=$m0, start_rep=$start_rep, B_REPS=$B_REPS"

Rscript run_simex_nhanes_serial_scaled.R $m0 $B_REPS $start_rep