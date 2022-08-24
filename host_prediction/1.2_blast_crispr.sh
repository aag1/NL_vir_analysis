#!/bin/bash
#SBATCH --job-name=job1.2
#SBATCH --output=job1.2_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/host_prediction



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -query /data/umcg-tifn/DATABASES/CRISPRCasdb/20210121_spacer_34_safe_names.fasta \
    -db DB3 \
    -dust no \
    -evalue 1 \
    -task blastn-short \
    -max_target_seqs 1000000 \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen' \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out blastn_spacers_CRISPRCasdb_vs_DB3.out



chmod 440 blastn_spacers_CRISPRCasdb_vs_DB3.out
