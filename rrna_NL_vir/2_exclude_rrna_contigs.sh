#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir

rm NL_vir_genome_fragments.n*



# ------------------------------ rrna contigs ids ------------------------------ #
awk '$8 - $7 + 1 > $14*0.5' NL_vir_vs_SILVA_lsu_rrna.txt | cut -d$'\t' -f2 | sort | uniq > NL_vir_lsu_contigs.ids

awk '$8 - $7 + 1 > $14*0.5' NL_vir_vs_SILVA_ssu_rrna.txt | cut -d$'\t' -f2 | sort | uniq > NL_vir_ssu_contigs.ids

cat NL_vir_lsu_contigs.ids > NL_vir_rrna_contigs.ids

comm -13 NL_vir_lsu_contigs.ids NL_vir_ssu_contigs.ids >> NL_vir_rrna_contigs.ids



# ------------------------------ exclude rrna contigs from FASTA ------------------------------ #
grep '>' /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta | sed 's/^>//' | sort > NL_vir_all_contigs.ids

comm -13 NL_vir_rrna_contigs.ids NL_vir_all_contigs.ids > NL_vir_sele_contigs.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta \
    NL_vir_sele_contigs.ids > NL_vir_genome_fragments.no_rrna.fasta



# ------------------------------ exclude rrna contigs from table ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/exclude_rrna_from_table.R



# ------------------------------ clean up & permissions ------------------------------ #
rm NL_vir_all_contigs.ids NL_vir_sele_contigs.ids

chmod 440 NL_vir_*_contigs.ids NL_vir_genome_fragments.no_rrna.*
