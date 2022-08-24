#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/DB4_in_GenBank



# https://www.ncbi.nlm.nih.gov/books/NBK179288/
# https://www.ncbi.nlm.nih.gov/books/NBK569846/

export PATH=${PATH}:${HOME}/edirect

/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -n Viruses
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -t 10239 > Viruses.txids

chmod 440 Viruses.txids
