#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=1gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info



module purge; module load R; module list

Rscript ${dir}/plot_TerL_msa.R



chmod 440 *png DB3_terL_msa_trimmed.pdf
