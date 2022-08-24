#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DEVoC_reads



module purge; module load R; module list

Rscript ${dir}/summarize_mapping_results_DB3.R



chmod 440 DEVoC_DB3_mapping_summary.txt
