#!/bin/bash
#SBATCH --job-name=job1.1
#SBATCH --output=job1.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /scratch/umcg-agulyaeva/DEVoC_reads



ARR=( $(awk -F'\t' '{if ($10 ~ /^HA/) print $7}' ${dir}/DEVoC_adult_cohort_PRJNA722819.tsv) )

for x in "${ARR[@]}"
do

    f1=$(echo ${x} | cut -d';' -f1)
    f2=$(echo ${x} | cut -d';' -f2)

    wget ${f1}
    wget ${f2}

done
