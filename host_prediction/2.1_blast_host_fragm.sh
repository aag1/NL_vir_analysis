#!/bin/bash
#SBATCH --job-name=job2.1
#SBATCH --output=job2.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/host_prediction




# ------------------------------ restrict blast taxa ------------------------------ #
# https://www.ncbi.nlm.nih.gov/books/NBK179288/
# https://www.ncbi.nlm.nih.gov/books/NBK569846/

export PATH=${PATH}:${HOME}/edirect

/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -n Procaryotae
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -t 2 > Procaryotae.txids
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -t 2157 >> Procaryotae.txids

chmod 440 Procaryotae.txids




# ------------------------------ list DB3 species with fragments cleaved off by Cenote-Taker2 ------------------------------ #
awk -F'\t' '{if ($7=="yes") print $1}' /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3_table1.txt | sort | uniq > DB3_vOTUs_with_prophage.txt

chmod 440 DB3_vOTUs_with_prophage.txt
