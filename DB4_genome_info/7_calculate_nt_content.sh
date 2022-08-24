#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_genome_info



/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta \
    DB4_plus.ids > DB4_plus.fasta



module purge; module load R; module list

Rscript ${dir}/calculate_nt_content.R \
    --seqF DB4_plus.fasta \
    --prefix DB4



chmod 440 DB4*
