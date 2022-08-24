#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot



module purge; module load R; module list

Rscript ${dir}/sele_prot_seq.R



module purge; module load MAFFT; module list

for x in 'RVT_1' 'Phage_integrase' 'rve'
do
    pfam_msa=$(ls /home/umcg-agulyaeva/NL_vir_analysis/DB4_proteome_annot/${x}_*_seed.fasta)

    mafft \
        --add 'DB4_'${x}'_seq.fasta' \
        --thread ${SLURM_CPUS_PER_TASK} \
        ${pfam_msa} > ${x}'_msa_with_pfam.fasta'
done



module purge; module load R; module list

Rscript ${dir}/trim_msa.R



rm *_seq.fasta *_msa_with_pfam.fasta

chmod 440 *_info.txt *_msa.fasta
