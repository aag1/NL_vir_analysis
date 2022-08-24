#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/benchmark_markers




# ------------------------------ explore group-specific VOG profiles ------------------------------ #
ARR1=( \
'Herviviricetes' 'Herpesviridae' 'Alloherpesviridae' 'Malacoherpesviridae' \
'Adenoviridae' \
'Papillomaviridae' \
'Polyomaviridae' \
'Preplasmiviricota' 'Maveriviricetes' 'Polintoviricetes' 'Tectiliviricetes' \
'Nucleocytoviricota' \
)


for g in "${ARR1[@]}"
do

    echo -e '\n'${g}

    ARR2=( $(grep ${g}'$' /data/umcg-tifn/DATABASES/VOG_207/vog.lca.tsv | cut -d$'\t' -f1) )

    for d in "${ARR2[@]}"
    do
        grep ${d} /data/umcg-tifn/DATABASES/VOG_207/vog.annotations.tsv
    done

    echo -e '\n'

done




# ------------------------------ explore VOG TerL profiles ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/explore_vog_terl_profiles.R

chmod 440 vog_terl_profiles.txt VOG00461_taxo.txt

echo -e '\n'
grep -v Caudoviricetes VOG00461_taxo.txt
