#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_in_GenBank




# ------------------------------ data for dotplots ------------------------------ #
module purge; module load EMBOSS; module list

awk -F'\t' '{if ($4 == "KX501134.1" || $9 ~ /Minot|Dzunkova/) print $1"\t"$4"\t"$11}' DB4_cognate_GenBank_genera_level_with_pub_info.txt > sele_pairs.txt

while read -r nl_id gb_id gb_rev
do
    echo ${nl_id} > ${nl_id}.id
    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.fasta ${nl_id}.id > ${nl_id}.fasta


    efetch -db nuccore -id ${gb_id} -mode text -format fasta > ${gb_id}.fasta
    if [ ${gb_rev} == 1 ]; then revseq ${gb_id}.fasta ${gb_id}.fasta; fi


    pair_id=${nl_id}'___'${gb_id}
    cat ${nl_id}.fasta > ${pair_id}.fasta
    cat ${gb_id}.fasta >> ${pair_id}.fasta


    polydot ${pair_id}.fasta -wordsize 12 -graph pdf -dumpfeat -outfeat ${pair_id}.gff
    mv polydot.pdf ${pair_id}.pdf

done < sele_pairs.txt

rm sele_pairs.txt *id *fasta




# ------------------------------ plot dotplots ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/dotplots.R

chmod 440 *gff *pdf
