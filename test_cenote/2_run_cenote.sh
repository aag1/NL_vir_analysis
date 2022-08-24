#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/test_cenote



N=$((${SLURM_ARRAY_TASK_ID}-1))
CHUNK_ID=viral_refseq_209_${N}



module purge; module load Anaconda3; module list
source activate cenote-taker2_env; conda list

python /data/umcg-agulyaeva/SOFTWARE/Cenote-Taker2/unlimited_breadsticks.py \
    --contigs ${CHUNK_ID}.fna \
    --run_title ${CHUNK_ID}_vir \
    --virus_domain_db 'virion' \
    --minimum_length_circular 3000 \
    --minimum_length_linear 10000 \
    --circ_minimum_hallmark_genes 1 \
    --lin_minimum_hallmark_genes 2 \
    --prune_prophage True \
    --filter_out_plasmids True \
    --mem 4 \
    --cpu ${SLURM_CPUS_PER_TASK}

conda deactivate



mv ${CHUNK_ID}_vir/${CHUNK_ID}_vir_CONTIG_SUMMARY.tsv .
mv ${CHUNK_ID}_vir/${CHUNK_ID}_vir_PRUNING_INFO_TABLE.tsv .

chmod 440 ${CHUNK_ID}_vir_CONTIG_SUMMARY.tsv
chmod 440 ${CHUNK_ID}_vir_PRUNING_INFO_TABLE.tsv

rm -rf ${CHUNK_ID}_vir
rm ${CHUNK_ID}.fasta
rm ${CHUNK_ID}.over_3000nt.fasta
