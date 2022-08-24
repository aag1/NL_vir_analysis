#!/bin/bash
#SBATCH --job-name=job1.0
#SBATCH --output=job1.0_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/host_prediction



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/makeblastdb \
    -dbtype 'nucl' \
    -in /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.fasta \
    -out DB3
