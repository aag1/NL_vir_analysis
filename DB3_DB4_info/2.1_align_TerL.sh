#!/bin/bash
#SBATCH --job-name=job2.1
#SBATCH --output=job2.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info




# -------------------------------------- list TerL ids -------------------------------------- #
module purge; module load R; module list

Rscript ${dir}/list_circ_caudo_terL.R




# -------------------------------------- get TerL from RefSeq -------------------------------------- #
awk -F'\t' '$3 == "RefSeq" && $2 != "" { print $2 }' DB2_circ_caudo_terL.txt > RefSeq_terL.ids


/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_protein.faa \
    RefSeq_terL.ids > DB2_circ_caudo_terL.fasta




# -------------------------------------- get TerL from NL_vir -------------------------------------- #
awk -F'\t' '$3 == "NL_vir" && $2 != "" { print $2 }' DB2_circ_caudo_terL.txt > NL_vir_terL.ids


/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.fasta \
    NL_vir_terL.ids >> DB2_circ_caudo_terL.fasta




# -------------------------------------- get TerL from other databases -------------------------------------- #
for db in 'MGV' 'GPD' 'GVD' 'DEVoC' 'HuVirDB' 'BanfieldLab'
do
    awk -v db=${db} -F'\t' '$3 == db && $2 != "" { print $2 }' DB2_circ_caudo_terL.txt > ${db}_terL.ids

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        /data/umcg-tifn/NL_vir_analysis/get_databases_caudo/${db}_5kb_cMAG_proteomes.fasta \
        ${db}_terL.ids >> DB2_circ_caudo_terL.fasta
done




rm *ids
chmod 440 DB2_circ_caudo_terL*
