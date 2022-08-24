#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/predict_proteomes



module purge; module load BBMap; module list

partition.sh \
    in=/data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_genome_fragments.no_rrna.fasta \
    out=NL_vir_genome_fragments_%_NT.fasta \
    ways=100
