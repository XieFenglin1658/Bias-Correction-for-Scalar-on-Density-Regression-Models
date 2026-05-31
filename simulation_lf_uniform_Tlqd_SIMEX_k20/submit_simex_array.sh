#!/bin/bash
#SBATCH --account=5sigma
#SBATCH --job-name=simex_multi
#SBATCH --output=logs/simex_%A_%a.out
#SBATCH --error=logs/simex_%A_%a.err
#SBATCH --time=0-10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-600%200

cd /insomnia001/depts/5sigma/users/fx2212/simulation_lf_Tlqd
mkdir -p logs

source /insomnia001/home/fx2212/miniconda3/bin/activate fda-43

n_list=(50 100 200)
sigma_list=(0.01 0.05)
rep_list=($(seq 1 100))

idx=$((SLURM_ARRAY_TASK_ID - 1))
n_sigma=${#sigma_list[@]}
n_rep=${#rep_list[@]}
total_per_n=$((n_sigma * n_rep))

n_index=$((idx / total_per_n))
remainder=$((idx % total_per_n))
sigma_index=$((remainder / n_rep))
rep_index=$((remainder % n_rep))

n=${n_list[$n_index]}
sigma=${sigma_list[$sigma_index]}
rep=$((rep_index + 1))

Rscript run_simex_single.R $n $sigma $rep