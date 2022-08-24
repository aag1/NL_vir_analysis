#!/bin/bash
#SBATCH --job-name=job4.2
#SBATCH --output=job4.2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-135
#SBATCH --export=NONE
#SBATCH --get-user-env=L



file=/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD_scripts_logs_info/LLD_raw_reads_number_sele.txt
N=$((${SLURM_ARRAY_TASK_ID} + 1000))
SAMPLE_ID=$(sed "${N}q;d" ${file} | cut -d$'\t' -f1)
echo '-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'



cd /data/umcg-tifn/NL_vir_analysis/run_cenote



module purge; module load Anaconda3; module list
source activate cenote-taker2_env; conda list

python /data/umcg-agulyaeva/SOFTWARE/Cenote-Taker2/unlimited_breadsticks.py \
    --contigs /data/umcg-tifn/LLD_assembly/${SAMPLE_ID}/${SAMPLE_ID}_contigs.fasta \
    --run_title ${SAMPLE_ID}_vir \
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



mv ${SAMPLE_ID}_vir/${SAMPLE_ID}_vir_CONTIG_SUMMARY.tsv .
mv ${SAMPLE_ID}_vir/${SAMPLE_ID}_vir_PRUNING_INFO_TABLE.tsv .
mv ${SAMPLE_ID}_vir/final_combined_virus_sequences_${SAMPLE_ID}_vir.fna .

chmod 440 ${SAMPLE_ID}_vir_CONTIG_SUMMARY.tsv
chmod 440 ${SAMPLE_ID}_vir_PRUNING_INFO_TABLE.tsv
chmod 440 final_combined_virus_sequences_${SAMPLE_ID}_vir.fna

rm -rf ${SAMPLE_ID}_vir
rm -rf ${SAMPLE_ID}_contigs.fasta
rm ${SAMPLE_ID}_contigs.over_3000nt.fasta
