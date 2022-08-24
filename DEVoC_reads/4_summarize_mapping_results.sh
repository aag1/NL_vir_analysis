#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DEVoC_reads



module purge; module load R; module list

Rscript ${dir}/summarize_mapping_results.R



chmod 440 DEVoC_DB1_mapping_summary.txt DEVoC_DB1_abundance.txt
