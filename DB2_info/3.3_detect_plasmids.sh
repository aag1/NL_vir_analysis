#!/bin/bash
#SBATCH --job-name=job3.3
#SBATCH --output=job3.3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB2_info



module purge; module load R; module list

Rscript ${dir}/collect_plasx_data.R

chmod 440 DB2_circ_PlasX_score.txt
