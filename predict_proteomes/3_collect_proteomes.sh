#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


cd /data/umcg-tifn/NL_vir_analysis/predict_proteomes


cat NL_vir_genome_fragments_*_proteomes.fasta > NL_vir_genome_fragments_proteomes.fasta


sed -n '1p' NL_vir_genome_fragments_0_proteomes.txt > NL_vir_genome_fragments_proteomes.txt
for i in {0..99}; do sed '1d' NL_vir_genome_fragments_${i}_proteomes.txt >> NL_vir_genome_fragments_proteomes.txt; done


sed -n '1p' NL_vir_genome_fragments_0_SeleCode.txt > NL_vir_genome_fragments_SeleCode.txt
for i in {0..99}; do sed '1d' NL_vir_genome_fragments_${i}_SeleCode.txt >> NL_vir_genome_fragments_SeleCode.txt; done


rm NL_vir_genome_fragments_*_*
chmod 440 *
