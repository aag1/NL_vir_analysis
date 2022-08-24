#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-54
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot




# ------------------------------ genome ID ------------------------------ #
GENOME_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids)

GENOME_SOURCE=$(awk -v var=${GENOME_ID} '{if ($1 == var) print $8}' /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt)

echo '#################### WORKING WITH '${GENOME_ID}' FROM '${GENOME_SOURCE}' ####################'




# ------------------------------ get proteome ------------------------------ #
if [ ${GENOME_SOURCE} == 'NL_vir' ]; then

    f1='/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt'
    f2='/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.fasta'

elif [ ${GENOME_SOURCE} == 'Benler_2020' ]; then

    f1='/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.txt'
    f2='/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.fasta'

else

    f1='/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/'${GENOME_SOURCE}'_5kb_cMAG_proteomes.txt'
    f2='/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/'${GENOME_SOURCE}'_5kb_cMAG_proteomes.fasta'

fi

awk -F'\t' -v var=${GENOME_ID} '{if ($1 == var) print $2}' ${f1} > ${GENOME_ID}_AA.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${f2} ${GENOME_ID}_AA.ids > ${GENOME_ID}_AA.fasta




# ------------------------------ run hmmer ------------------------------ #
module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --domtblout ${GENOME_ID}_vs_Pfam.txt \
    -o ${GENOME_ID}_vs_Pfam.out \
    /data/umcg-tifn/DATABASES/PfamA_35.0/Pfam-A.hmm \
    ${GENOME_ID}_AA.fasta




# ------------------------------ parse hmmer ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/parse_hmmer.R ${GENOME_ID}




rm ${GENOME_ID}_AA* ${GENOME_ID}_vs_Pfam.out

chmod 440 ${GENOME_ID}_vs_Pfam.txt ${GENOME_ID}_annot.txt
