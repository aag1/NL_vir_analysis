#!/bin/bash
#SBATCH --job-name=job2.1
#SBATCH --output=job2.1_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=12:00:00
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



### fix fastq read identifiers (see https://forum.biobakery.org/t/kneaddata-fastq-header-problem/482)
for i in 1 2
do
    gzip -cd ${SAMPLE_ID}_${i}.fastq.gz | sed 's/^\(@[A-Z0-9]\+\.[0-9]\+\) [0-9]\+\(\/[12]\)$/\1\2/' | gzip > ${SAMPLE_ID}_${i}.fixed.fastq.gz
done



### raw reads quality
module purge; module load FastQC; module list

fastqc \
    --threads 2 \
    ${SAMPLE_ID}_1.fixed.fastq.gz \
    ${SAMPLE_ID}_2.fixed.fastq.gz



### clean reads by KneadData
module purge; module load Anaconda3; module list

source activate KneadData
conda list

trimmomatic_path=/home/umcg-agulyaeva/SOFTWARE/Trimmomatic-0.33

kneaddata \
    --input ${SAMPLE_ID}_1.fixed.fastq.gz \
    --input ${SAMPLE_ID}_2.fixed.fastq.gz \
    --reference-db /data/umcg-tifn/bowtie2_human_genome_reference \
    --output-prefix ${SAMPLE_ID}_kneaddata \
    --output ${SAMPLE_ID}_kneaddata_out \
    --trimmomatic ${trimmomatic_path} \
    --trimmomatic-options "ILLUMINACLIP:${trimmomatic_path}/adapters/NexteraPE-PE.fa:2:30:10:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --threads ${SLURM_CPUS_PER_TASK}

conda deactivate



### remove extra files
mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata_paired_1.fastq .
mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata_paired_2.fastq .
mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata.log .
rm -rf ${SAMPLE_ID}_kneaddata_out



### clean reads quality
module purge; module load FastQC; module list

fastqc \
    --threads 2 \
    ${SAMPLE_ID}_kneaddata_paired_1.fastq \
    ${SAMPLE_ID}_kneaddata_paired_2.fastq



### compress clean reads files
gzip ${SAMPLE_ID}_kneaddata_paired_1.fastq
gzip ${SAMPLE_ID}_kneaddata_paired_2.fastq
