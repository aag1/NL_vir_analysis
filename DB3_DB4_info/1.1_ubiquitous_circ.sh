#!/bin/bash
#SBATCH --job-name=job1.1
#SBATCH --output=job1.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info




# ------------------------------ get sequence ------------------------------ #
sed '1d' /data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt | awk -F'\t' '$9 == "TR" && $2 > 500' | cut -d$'\t' -f1 > ubiquitous_circ.id

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta \
    ubiquitous_circ.id > ubiquitous_circ.fasta

rm ubiquitous_circ.id
chmod 440 ubiquitous_circ.fasta




# ------------------------------ run blast ------------------------------ #
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -task 'blastn' \
    -db 'nt' \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query ubiquitous_circ.fasta \
    -out ubiquitous_circ_vs_NT.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident staxids sscinames'

chmod 440 ubiquitous_circ_vs_NT.txt

awk '$4 > 5000' ubiquitous_circ_vs_NT.txt
