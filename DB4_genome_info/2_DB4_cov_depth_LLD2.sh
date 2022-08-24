#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-338
#SBATCH --export=NONE
#SBATCH --get-user-env=L



file=/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD2_scripts_logs_info/LLD2_raw_reads_number.txt
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file} | cut -d$'\t' -f1)
echo '#################### WORKING WITH '${SAMPLE_ID}' SAMPLE ####################'



cd /data/umcg-tifn/NL_vir_analysis/DB4_genome_info



module purge; module load SAMtools; module list

samtools depth /data/umcg-tifn/NL_vir_analysis/map_reads/map_LLD2_reads/${SAMPLE_ID}/${SAMPLE_ID}.sorted.bam |
grep -f DB4_plus.ids > ${SAMPLE_ID}_DB4_cov_depth.txt

chmod 440 ${SAMPLE_ID}_DB4_cov_depth.txt