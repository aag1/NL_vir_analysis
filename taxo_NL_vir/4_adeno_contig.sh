#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/taxo_NL_vir



# ------------------------------ get sequence ------------------------------ #
awk -F'\t' '{if ($11 == "Adenoviridae") print $9}' NL_vir_genome_fragments.no_rrna.taxo.txt > NL_vir_adeno.id

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta \
    NL_vir_adeno.id > NL_vir_adeno.fasta

rm NL_vir_adeno.id
chmod 440 NL_vir_adeno.fasta



# ------------------------------ run blast ------------------------------ #
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -task 'blastn' \
    -db 'nt' \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query NL_vir_adeno.fasta \
    -out NL_vir_adeno_vs_NT.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident staxids sscinames'

chmod 440 NL_vir_adeno_vs_NT.txt
