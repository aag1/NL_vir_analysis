#!/bin/bash
#SBATCH --job-name=job3.1
#SBATCH --output=job3.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/DB2_info



# ------------------------------ extract vOTU representatives with terminal repeats from DB2 ------------------------------ #
sed '1d' DB2_info.txt | awk -F'\t' '{if ($9 == "TR") print $1}' > DB2_circ.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1.fasta \
    DB2_circ.ids > DB2_circ.fasta

rm DB2_circ.ids



# ------------------------------ split multifasta into 10 parts ------------------------------ #
module purge; module load BBMap; module list

partition.sh \
    in=DB2_circ.fasta \
    out=DB2_circ_%.fasta \
    ways=10

chmod 440 DB2_circ*.fasta
