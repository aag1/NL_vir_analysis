#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/benchmark_markers



module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --noali \
    --domtblout markers_vs_viral_refseq.txt \
    -o hmmer_out.txt \
    --cpu ${SLURM_CPUS_PER_TASK} \
    markers.hmm \
    /data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_protein.faa



rm hmmer_out.txt

chmod 440 markers_vs_viral_refseq.txt
