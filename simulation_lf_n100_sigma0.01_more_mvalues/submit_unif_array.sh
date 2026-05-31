#!/bin/bash
#SBATCH --account=5sigma
#SBATCH --job-name=unif_lf
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=0-08:00          
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-12             
#SBATCH --output=unif_%a_%j.out
#SBATCH --error=unif_%a_%j.err

source /insomnia001/home/fx2212/miniconda3/bin/activate fda-43
cd /insomnia001/depts/5sigma/users/fx2212/simulation_lf_n100_sigma0.01_more_mvalues


M_LIST=(50 100 200 300 500 1000 2000 3000 5000 10000 100000 1000000)
M_TARGET=${M_LIST[$SLURM_ARRAY_TASK_ID-1]}

Rscript run_unif_lf.R $M_TARGET 100