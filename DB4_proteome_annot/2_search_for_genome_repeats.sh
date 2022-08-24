#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot



while read genome_id
do

    echo ${genome_id} > ${genome_id}.id


    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.fasta \
        ${genome_id}.id > ${genome_id}.fa


    /home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
        -task 'blastn' \
        -query ${genome_id}.fa \
        -subject ${genome_id}.fa \
        -out ${genome_id}_repeats.11.txt \
        -outfmt 11 \
        -evalue 0.001


    /home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blast_formatter \
        -archive ${genome_id}_repeats.11.txt \
        -out ${genome_id}_repeats.0.txt \
        -outfmt 0


    /home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blast_formatter \
        -archive ${genome_id}_repeats.11.txt \
        -out ${genome_id}_repeats.6.txt \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident'

done < /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids



module purge; module load R; module list

Rscript ${dir}/summarize_repeats.R



rm *.id *.fa
chmod 440 *_repeats*txt
