#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info




# -------------------------------------- align TerL MSAs -------------------------------------- #
module purge; module load GCCcore/6.4.0; module list


/home/umcg-agulyaeva/SOFTWARE/hh-suite/build/bin/hhalign \
    -M 50 \
    -mact 0 \
    -all \
    -i /data/umcg-tifn/DATABASES/VOG_207/vog_msa/VOG00461.msa \
    -t /data/umcg-tifn/DATABASES/data_Benler_2020/alignments_and_newick_files/Fig1_TerL.fasta \
    -o TerL_2x.hhr \
    -oa3m TerL_2x.a3m

/home/umcg-agulyaeva/SOFTWARE/hh-suite/build/scripts/reformat.pl \
    TerL_2x.a3m \
    TerL_2x.fasta \
    -uc \
    -g '-' \
    -l 80


/home/umcg-agulyaeva/SOFTWARE/hh-suite/build/bin/hhalign \
    -M 50 \
    -mact 0 \
    -all \
    -i TerL_2x.fasta \
    -t /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/trees/TerL_cMAG.afa \
    -o TerL_3x.hhr \
    -oa3m TerL_3x.a3m

/home/umcg-agulyaeva/SOFTWARE/hh-suite/build/scripts/reformat.pl \
    TerL_3x.a3m \
    TerL_3x.fasta \
    -uc \
    -g '-' \
    -l 80




# -------------------------------------- add sequences -------------------------------------- #
module purge; module load MAFFT; module list
mafft \
    --add DB2_circ_caudo_terL.fasta \
    --thread ${SLURM_CPUS_PER_TASK} \
    TerL_3x.fasta > TerL_msa.fasta


rm TerL_2x* TerL_3x.hhr TerL_3x.a3m
chmod 440 TerL_3x.fasta TerL_msa.fasta
