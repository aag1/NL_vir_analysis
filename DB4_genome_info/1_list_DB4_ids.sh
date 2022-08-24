#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


cd /data/umcg-tifn/NL_vir_analysis/DB4_genome_info


sed '1d' /data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt |
awk -F'\t' '{if (($9 == "TR") && ($8 == "Caudoviricetes") && (($2 / 1135 > 0.05) || ($3 / 338 > 0.05) || ($4 / 298 > 0.05) || ($5 / 455 > 0.05))) print $1}' > DB4_plus.ids


chmod 440 DB4_plus.ids
