#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-395
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/host_prediction




# ------------------------------ species representative ID ------------------------------ #
sp_id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" DB3_vOTUs_with_prophage.txt)
echo '#################### WORKING WITH '${sp_id}' CLUSTER ####################'




# ------------------------------ get host fragments cleaved off by Cenote-Taker2 ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/get_host_fragm.R ${sp_id}




# ------------------------------ blast host fragments against nr procaryotes ------------------------------ #
/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -task 'blastn' \
    -db 'nt' \
    -taxidlist 'Procaryotae.txids' \
    -perc_identity 95 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -query ${sp_id}'_host_fragm.fasta' \
    -out ${sp_id}'_host_fragm_vs_ncbi.txt' \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident staxids sscinames'

chmod 440 ${sp_id}_host_fragm*
