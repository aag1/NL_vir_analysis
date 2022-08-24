#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ extract potential NCLVD genomes ------------------------------ #
cd /data/umcg-tifn/NL_vir_analysis/taxo_NL_vir


awk '{if ($4 ~ /^GVOGm/) print $1}' markers_vs_NL_vir.txt | sed 's/_[0-9]\+_c[0-9]\+$//' | sort | uniq > potential_NCLVD_genomes.ids


/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta \
    potential_NCLVD_genomes.ids > potential_NCLVD_genomes.fasta




# ------------------------------ run ViralRecall ------------------------------ #
module purge
module load Python/3.9.5-GCCcore-10.3.0 Biopython/1.79-foss-2021a matplotlib/3.4.2-foss-2021a
module load prodigal/2.6.3-GCCcore-10.3.0
module load HMMER/3.3.2-gompi-2021a
module list


cp -pr /data/umcg-agulyaeva/SOFTWARE/viralrecall/hmm .
cp -pr /data/umcg-agulyaeva/SOFTWARE/viralrecall/acc .


python /data/umcg-agulyaeva/SOFTWARE/viralrecall/viralrecall.py \
    --input potential_NCLVD_genomes.fasta \
    --project nclvd_ViralRecall \
    --contiglevel \
    --evalue 1e-3 \
    --cpus ${SLURM_CPUS_PER_TASK}


rm potential_NCLVD_genomes.fasta out.txt err.txt
rm -rf hmm acc
chmod 750 nclvd_ViralRecall
chmod 440 potential_NCLVD_genomes.ids nclvd_ViralRecall/*
