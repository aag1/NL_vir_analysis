#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/test_cenote



module purge; module load BBMap; module list

partition.sh \
    in=/data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_genomic.fna \
    out=viral_refseq_209_%.fna \
    ways=100
