#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-54
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/DB4_in_BanfieldLab




# ------------------------------ genome ID ------------------------------ #
id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids)
echo '#################### WORKING WITH '${id}' CLUSTER ####################'




# ------------------------------ get genome ------------------------------ #
echo ${id} > ${id}.id

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.fasta \
    ${id}.id > ${id}.fasta

rm ${id}.id




# ------------------------------ run blastn ------------------------------ #
for db in 'Devoto_2019' 'Al-Shayeb_2020' 'Borges_2021'
do
    /home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
        -task 'blastn' \
        -db ${db} \
        -evalue 0.001 \
        -perc_identity 50 \
        -num_threads ${SLURM_CPUS_PER_TASK} \
        -query ${id}.fasta \
        -out ${id}_vs_${db}.txt \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident'
done

chmod 440 ${id}_vs_*.txt
rm ${id}.fasta
