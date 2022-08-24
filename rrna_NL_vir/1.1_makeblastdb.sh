#!/bin/bash
#SBATCH --job-name=job1.1
#SBATCH --output=job1.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir



module purge; module load BLAST+; module list

makeblastdb \
    -dbtype 'nucl' \
    -in /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta \
    -out NL_vir_genome_fragments
