#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/DEVoC_reads



module purge; module load Anaconda3; module list

source activate MultiQC
conda list

multiqc \
    --module fastqc \
    --filename DEVoC_raw_reads_multiqc_report.html \
    /scratch/umcg-agulyaeva/DEVoC_reads/*fixed_fastqc.zip

multiqc \
    --module fastqc \
    --filename DEVoC_clean_reads_multiqc_report.html \
    /scratch/umcg-agulyaeva/DEVoC_reads/*_kneaddata_paired_*_fastqc.zip

conda deactivate



chmod 440 DEVoC_*_reads_multiqc_report.html
