#!/bin/bash
#SBATCH --job-name=job1.3
#SBATCH --output=job1.3_%A.out
#SBATCH --mem=20gb
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir



module purge; module load BLAST+; module list

blastn \
    -task 'blastn' \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -db NL_vir_genome_fragments \
    -query /data/umcg-tifn/DATABASES/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
    -out 'NL_vir_vs_SILVA_ssu_rrna.txt' \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen'



chmod 440 NL_vir_vs_SILVA_ssu_rrna.txt
