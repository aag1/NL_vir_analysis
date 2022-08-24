#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/taxo_NL_vir



module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --noali \
    --domtblout markers_vs_NL_vir.txt \
    -o hmmer_out.txt \
    --cpu ${SLURM_CPUS_PER_TASK} \
    /data/umcg-tifn/NL_vir_analysis/benchmark_markers/markers.hmm \
    /data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.fasta



rm hmmer_out.txt
chmod 440 markers_vs_NL_vir.txt
