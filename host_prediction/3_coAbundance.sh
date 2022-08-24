#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/host_prediction



module purge; module load R; module list

Rscript ${dir}/coAbundance.R

Rscript ${dir}/compare_host_predictions.R



chmod 440 DB4_phage_host_coAbundance_metaAnalysis* DB4_compare_host_predictions.txt
