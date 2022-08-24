#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/test_code_prediction




sed '1d' /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt | cut -d$'\t' -f1 | sort | uniq > crAss378.txt

while read GENOME_ID
do

    echo -e '\n\n\nWorking with '${GENOME_ID}' ...'


    # get genome
    echo ${GENOME_ID} > ${GENOME_ID}.id

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.fasta \
        ${GENOME_ID}.id > ${GENOME_ID}.fasta


    # scan for tRNA genes
    /home/umcg-agulyaeva/SOFTWARE/trnascan_configure_dir/bin/tRNAscan-SE \
        -B \
        --output ${GENOME_ID}_tRNA_genes.txt \
        --thread ${SLURM_CPUS_PER_TASK} \
        --forceow \
        ${GENOME_ID}.fasta

done < crAss378.txt




module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/DB3_tRNA_scan/summarize_tRNAscan.R crAss378_tRNAscan.txt




rm crAss378.txt *id *fasta *_tRNA_genes.txt

chmod 440 crAss378_tRNAscan.txt
