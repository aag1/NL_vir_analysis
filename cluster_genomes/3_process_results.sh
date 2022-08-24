#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/cluster_genomes




# ------------------------------ genome lengths & GC content ------------------------------ #
module purge; module load EMBOSS; module list

infoseq \
    -auto \
    -nocolumns -delimiter $'\t' \
    -only -name -length -pgc \
    DB0.fasta > DB0_lengths_GC.txt

chmod 440 DB0_lengths_GC.txt




# ------------------------------ Benler et al. TerL IDs ------------------------------ #
grep '>' /data/umcg-tifn/DATABASES/data_Benler_2020/alignments_and_newick_files/Fig1_TerL.fasta | sed 's/^>//' | sort | uniq > all_TerL.ids

cut -d',' -f2 /data/umcg-tifn/DATABASES/data_Benler_2020/gut_phage_annotations.csv | sort | uniq > gut_proteome.ids

comm -12 all_TerL.ids gut_proteome.ids | sed 's/^\([^_]\+\)_\([0-9]\+\)$/\1\t\1_\2/' > Benler_2020_on_Fig1_TerL.txt

rm all_TerL.ids gut_proteome.ids

chmod 440 Benler_2020_on_Fig1_TerL.txt




# ------------------------------ BanfieldLab phage IDs by paper ------------------------------ #
grep '>' /data/umcg-tifn/DATABASES/data_Devoto_2019/15_Lak_phage_genomes.fasta |
sed 's/^>//' > 15_Lak_phage_genomes.ids

grep '>' /data/umcg-tifn/DATABASES/data_Al-Shayeb_2020/349_huge_phage_genomes.fasta |
sed 's/^>//' > 349_huge_phage_genomes.ids

chmod 440 15_Lak_phage_genomes.ids 349_huge_phage_genomes.ids




# ------------------------------ summarize clustering results ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/process_results.R


sed '1d' DB1_info.txt | cut -d$'\t' -f16 | sort | uniq > DB1.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    DB0.fasta \
    DB1.ids > DB1.fasta


rm DB1.ids

chmod 440 DB1_info.txt DB1.fasta
