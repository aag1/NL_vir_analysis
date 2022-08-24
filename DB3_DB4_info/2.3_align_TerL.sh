#!/bin/bash
#SBATCH --job-name=job2.3
#SBATCH --output=job2.3_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info



module purge; module load R; module list

Rscript ${dir}/trim_msa.R



chmod 440 DB3_terL_msa_*.fasta DB3.ids
