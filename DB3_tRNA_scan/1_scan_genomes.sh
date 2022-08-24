#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_tRNA_scan




while read GENOME_ID
do

    echo -e '\n\n\nWorking with '${GENOME_ID}' ...'


    # get genome
    echo ${GENOME_ID} > ${GENOME_ID}.id

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.fasta \
        ${GENOME_ID}.id > ${GENOME_ID}.fasta


    # scan for tRNA genes
    /home/umcg-agulyaeva/SOFTWARE/trnascan_configure_dir/bin/tRNAscan-SE \
        -B \
        --output ${GENOME_ID}_tRNA_genes.txt \
        --thread ${SLURM_CPUS_PER_TASK} \
        --forceow \
        ${GENOME_ID}.fasta

done < /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids




module purge; module load R; module list

Rscript ${dir}/summarize_tRNAscan.R DB3_tRNAscan.txt




rm *id *fasta *_tRNA_genes.txt

chmod 440 DB3_tRNAscan.txt
