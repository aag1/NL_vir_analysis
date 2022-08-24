#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info



module purge; module load R; module list

Rscript ${dir}/define_db3_db4.R



/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta \
    DB3.tree.ids > DB3.fasta



/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta \
    DB4.tree.ids > DB4.fasta



/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_genome_fragments.no_rrna.fasta \
    NL_vir_from_DB3_vOTUs.ids > NL_vir_from_DB3_vOTUs.fasta



rm NL_vir_from_DB3_vOTUs.ids

chmod 440 *.tree.ids DB3.fasta DB4.fasta DB3_table*.txt NL_vir_from_DB3_vOTUs.fasta
