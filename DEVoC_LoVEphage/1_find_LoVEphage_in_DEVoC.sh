#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ get LoVEphage genomes ------------------------------ #
if [ -f LoVEphage.fasta ]
then
    rm -rf LoVEphage.fasta
fi


for gb_id in 'MW660583' 'MZ919976' 'MZ919981' 'MZ919987'
do
    efetch -db nuccore -id ${gb_id} -mode text -format fasta >> LoVEphage.fasta
done




# ------------------------------ makeblastdb DEVoC ------------------------------ #
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/makeblastdb \
    -in /data/umcg-tifn/DATABASES/DEVoC/DEVoC.fasta \
    -dbtype nucl \
    -out DEVoC




# ------------------------------ run blast ------------------------------ #
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -task 'blastn' \
    -db DEVoC \
    -perc_identity 95 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query LoVEphage.fasta \
    -out LoVEphage_vs_DEVoC.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident'


echo -e '\n\n\n'
awk -F'\t' '$4 > 10000' LoVEphage_vs_DEVoC.txt


echo -e '\n\n\n'
grep ViralSequence127 /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3_table1.txt


rm DEVoC*
chmod 440 LoVEphage.fasta LoVEphage_vs_DEVoC.txt
