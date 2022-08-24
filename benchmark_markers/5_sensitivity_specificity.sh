#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/benchmark_markers



module purge; module load R; module list

Rscript ${dir}/select_TerL_hits.R \
    --hmmerF markers_vs_viral_refseq.txt \
    --prefix viral_refseq

Rscript ${dir}/sensitivity_specificity.R



chmod 440 viral_refseq_TerL_*.txt *_FP.txt *_FN.txt
