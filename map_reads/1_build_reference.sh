#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/map_reads

file=/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta



### build bowtie index
module purge; module load Bowtie2; module list

bowtie2-build \
    ${file} \
    DB1 \
    --threads ${SLURM_CPUS_PER_TASK}



### generate BED file
module purge; module load EMBOSS; module list

infoseq \
    -auto \
    -nocolumns -noheading \
    -delimiter $'\t0\t' \
    -only -name -length \
    ${file} > DB1.bed

chmod 440 DB1.bed
