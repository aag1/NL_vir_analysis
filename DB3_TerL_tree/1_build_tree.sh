#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_TerL_tree



module purge; module load Anaconda3; module list
source activate IQ-TREE

iqtree \
    -B 1000 \
    -mset WAG \
    -s /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3_terL_msa_trimmed.fasta \
    --prefix DB3_terL_tree \
    -T ${SLURM_CPUS_PER_TASK}

conda deactivate

chmod 440 DB3_terL_tree*
