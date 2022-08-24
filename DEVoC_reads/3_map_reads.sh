#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A_%a.out
#SBATCH --mem=50gb
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-254
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /scratch/umcg-agulyaeva/DEVoC_reads



### sample ID
if [ "${SLURM_ARRAY_TASK_ID}" -le 104 ]
then
    N=$((${SLURM_ARRAY_TASK_ID} + 1))
    file1=$(sed "${N}q;d" ${dir}/DEVoC_adult_cohort_PRJNA722819.tsv | cut -d$'\t' -f7 | cut -d';' -f1)
else
    N=$((${SLURM_ARRAY_TASK_ID} - 104 + 1))
    file1=$(sed "${N}q;d" ${dir}/DEVoC_pediatric_cohort_PRJNA723467.tsv | cut -d$'\t' -f7 | cut -d';' -f1)
fi

SAMPLE_ID=$(basename ${file1} _1.fastq.gz)

echo '#################### WORKING WITH '${SAMPLE_ID}' SAMPLE ####################'



### align reads
dir1=/data/umcg-tifn/NL_vir_analysis/map_reads

dir2=/data/umcg-tifn/NL_vir_analysis/DEVoC_reads

module purge; module load Bowtie2 SAMtools BEDTools; module list

bowtie2 \
    --very-sensitive \
    -x ${dir1}/DB1 \
    -1 ${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${SAMPLE_ID}.sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${SAMPLE_ID}.sorted.bam



### calculate statistics
bedtools coverage -a ${dir1}/DB1.bed -b ${SAMPLE_ID}.sorted.bam > ${dir2}/${SAMPLE_ID}.coverage.txt

chmod 440 ${dir2}/${SAMPLE_ID}.coverage.txt
