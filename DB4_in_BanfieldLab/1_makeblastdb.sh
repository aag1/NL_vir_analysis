#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/DB4_in_BanfieldLab



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/makeblastdb \
    -in /data/umcg-tifn/DATABASES/data_Devoto_2019/15_Lak_phage_genomes.fasta \
    -dbtype nucl \
    -out Devoto_2019



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/makeblastdb \
    -in /data/umcg-tifn/DATABASES/data_Al-Shayeb_2020/349_huge_phage_genomes.fasta \
    -dbtype nucl \
    -out Al-Shayeb_2020



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/makeblastdb \
    -in /data/umcg-tifn/DATABASES/data_Borges_2021_v1/borges_et_al_2021_phage_genomes.fasta \
    -dbtype nucl \
    -out Borges_2021
