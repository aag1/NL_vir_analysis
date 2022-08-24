#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/benchmark_markers

hmm_file='markers.hmm'
if [ -f ${hmm_file} ]; then rm -rf ${hmm_file}; fi




# ------------------------------ build TerL profiles ------------------------------ #
module purge; module load HMMER; module list

hmmbuild -n 'TerL_caudo' TerL_caudo.hmm /data/umcg-tifn/DATABASES/data_Benler_2020/alignments_and_newick_files/Fig1_TerL.fasta

hmmbuild -n 'TerL_crAss' TerL_crAss.hmm /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/trees/TerL_all.afa

hmmbuild -n 'TerL_VOG00461' TerL_VOG00461.hmm /data/umcg-tifn/DATABASES/VOG_207/vog_msa/VOG00461.msa


cat TerL_caudo.hmm > ${hmm_file}
cat TerL_crAss.hmm >> ${hmm_file}
cat TerL_VOG00461.hmm >> ${hmm_file}


chmod 440 TerL*hmm




# ------------------------------ get MSAs for other profiles ------------------------------ #
# Already in ${dir}/marker_MSAs/ :
# Herpes_MCP_PF03122_full.fasta, Polyoma_coat_PF00718_full.fasta (from Pfam 35.0)
# GVOGm*aln (from https://faylward.github.io/GVDB/)


# Adenoviridae
cp -p \
    /data/umcg-tifn/DATABASES/VOG_207/vog_msa/VOG05391.msa \
    ${dir}/marker_MSAs/Adeno_Hexon_protein_VOG05391.fa


# Papillomaviridae
cp -p \
    /data/umcg-tifn/DATABASES/VOG_207/vog_msa/VOG05075.msa \
    ${dir}/marker_MSAs/Papilloma_MCP_L1_VOG05075.fa


# Tectiliviricetes
for x in 'Bam_Toil_MCP' 'FLiP_group_MCP' 'Gemmatimonas_MCP' 'PM2_MCP' 'PRD1_MCP' 'STIV_MCP' 'odin_MCP'
do
    wget \
        ftp://ftp.ncbi.nlm.nih.gov/pub/yutinn/DJR_MCP_2017/MCP_aln_trees/${x}.afa \
        -P ${dir}/marker_MSAs/
done


chmod 440 ${dir}/marker_MSAs/*.fa
chmod 440 ${dir}/marker_MSAs/*.afa




# ------------------------------ build other profiles ------------------------------ #
for f in ${dir}/marker_MSAs/*
do

    b=$(basename ${f} | sed -E 's/\.[a-z]+$//')

    if [ "${b}" == 'GVOGm0890_A32' ] || [ "${b}" == 'GVOGm0013_SFII' ]; then continue; fi

    hmmbuild -n ${b} ${b}.hmm ${f}

    cat ${b}.hmm >> ${hmm_file}

    chmod 440 ${b}.hmm

done


chmod 440 ${hmm_file}
