#!/bin/bash
#SBATCH --job-name=job3.2
#SBATCH --output=job3.2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-10
#SBATCH --export=NONE
#SBATCH --get-user-env=L




N=$((${SLURM_ARRAY_TASK_ID}-1))

PREFIX='DB2_circ_'${N}

THREADS=${SLURM_CPUS_PER_TASK}

# see https://github.com/michaelkyu/PlasX




# ------------------------------ Step 1. Identify genes and annotate COGs and Pfams ------------------------------ #
module purge; module load Anaconda3; module list
source activate anvio-7.1


# remove additional info from name lines (e.g. virus name for RefSeq entries)
sed 's/^\(>[^ ]\+\) .\+$/\1/' $PREFIX.fasta > $PREFIX.simple_names.fasta


# Create an anvio contigs database from the fasta file
# - The `-L 0` parameter ensures that contigs remain intact and aren't split
anvi-gen-contigs-database \
    -L 0 \
    -T $THREADS \
    --project-name $PREFIX \
    -f $PREFIX.simple_names.fasta \
    -o $PREFIX.db


# Export gene calls (including amino acid sequences) to text file
anvi-export-gene-calls \
    --gene-caller prodigal \
    -c $PREFIX.db \
    -o $PREFIX-gene-calls.txt


# Annotate COGs
anvi-run-ncbi-cogs \
    -T $THREADS \
    --cog-version COG14 \
    --cog-data-dir /groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/anvio_COG_2014 \
    -c $PREFIX.db


# Annotate Pfams
anvi-run-pfams \
    -T $THREADS \
    --pfam-data-dir /groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/anvio_Pfam_v32 \
    -c $PREFIX.db


# Export functions to text file
anvi-export-functions \
    --annotation-sources COG14_FUNCTION,Pfam \
    -c $PREFIX.db \
    -o $PREFIX-cogs-and-pfams.txt


conda deactivate




# ------------------------------ Step 2. Annotate de novo gene families ------------------------------ #
source activate plasx

plasx search_de_novo_families \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-de-novo-families.txt \
    --threads $THREADS \
    --splits 32 \
    --overwrite




# ------------------------------ Step 3. Classify contigs as plasmid or non-plasmid ------------------------------ #
plasx predict \
    -a $PREFIX-cogs-and-pfams.txt $PREFIX-de-novo-families.txt \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-scores.txt \
    --overwrite

conda deactivate

chmod 440 $PREFIX*
