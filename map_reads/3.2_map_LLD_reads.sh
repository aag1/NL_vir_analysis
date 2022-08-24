#!/bin/bash
#SBATCH --job-name=job3.2
#SBATCH --output=job3.2_%A_%a.out
#SBATCH --mem=10gb
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-135
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### sample ID
file=/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD_scripts_logs_info/LLD_raw_reads_number_sele.txt
N=$((${SLURM_ARRAY_TASK_ID} + 1000))
SAMPLE_ID=$(sed "${N}q;d" ${file} | cut -d$'\t' -f1)
echo '#################### WORKING WITH '${SAMPLE_ID}' SAMPLE ####################'



### output directory
dir1=/data/umcg-tifn/NL_vir_analysis/map_reads

dir2=/data/umcg-tifn/NL_vir_analysis/map_reads/map_LLD_reads/${SAMPLE_ID}

if [ -d ${dir2} ]; then rm -rf ${dir2}; fi
mkdir ${dir2}
chmod 750 ${dir2}
cd ${dir2}



### align reads
module purge; module load Bowtie2 SAMtools BEDTools; module list

bowtie2 \
    --very-sensitive \
    -x ${dir1}/DB1 \
    -1 /data/umcg-tifn/LLD_clean_reads/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 /data/umcg-tifn/LLD_clean_reads/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${SAMPLE_ID}.sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${SAMPLE_ID}.sorted.bam



### calculate statistics
bedtools coverage -a ${dir1}/DB1.bed -b ${SAMPLE_ID}.sorted.bam > ${SAMPLE_ID}.coverage.txt

chmod 440 *
