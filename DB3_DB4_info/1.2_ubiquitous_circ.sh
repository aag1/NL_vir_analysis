#!/bin/bash
#SBATCH --job-name=job1.2
#SBATCH --output=job1.2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info



module purge; module load R; module list

Rscript ${dir}/plot_ubiquitous_circ.R



chmod 440 NL_vir005341_genome_map.pdf
