#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=40gb
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/cluster_genomes




# ------------------------------ Viral RefSeq ------------------------------ #
db='RefSeq'



### selected genomes
sed '1d' /data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_taxo.txt |
grep -v Riboviria |
awk -F'\t' '$3 > 2999' |
cut -d$'\t' -f1 > ${db}.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_genomic.fna \
    ${db}.ids > ${db}.fasta

chmod 440 ${db}.ids



### genomes with terminal repeats
perl /home/umcg-agulyaeva/NL_vir_analysis/taxo_NL_vir/identify_term_rep.pl \
    ${db}.fasta \
    ${db}_cMAG.txt \
    20

chmod 440 ${db}_cMAG.txt



### genetic code
module purge; module load prodigal; module list

for n in 11 4 15
do
    prodigal \
        -g ${n} \
        -p meta \
        -f gff \
        -o ${db}_AA_c${n}.gff \
        -a ${db}_AA_c${n}.fasta \
        -i ${db}.fasta

    sed -i "s/^\(>[^ ]\+\)/\1_c${n}/" ${db}_AA_c${n}.fasta
done


module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/predict_proteomes/select_code.R --prefix ${db}_AA


rm ${db}_AA_c*

chmod 440 ${db}_AA*



### add selected genomes to DB0
cat ${db}.fasta > DB0.fasta

rm ${db}.fasta




# ------------------------------ Benler et al. 2020 ------------------------------ #
db='Benler_2020'



### selected genomes
grep 'Uroviricota' /data/umcg-tifn/DATABASES/data_Benler_2020/40168_2021_1017_MOESM3_ESM.csv | cut -d',' -f1 > ${db}.ids

readarray arr < ${db}.ids

for g in ${arr[@]}
do
    f=/data/umcg-tifn/DATABASES/data_Benler_2020/nucl/${g}.fna
    cat ${f} >> ${db}.fasta
done

chmod 440 ${db}.ids



### genomes with terminal repeats
perl /home/umcg-agulyaeva/NL_vir_analysis/taxo_NL_vir/identify_term_rep.pl \
    ${db}.fasta \
    ${db}_cMAG.txt \
    20

chmod 440 ${db}_cMAG.txt



### genetic code
module purge; module load prodigal; module list

for n in 11 4 15
do
    prodigal \
        -g ${n} \
        -p meta \
        -f gff \
        -o ${db}_AA_c${n}.gff \
        -a ${db}_AA_c${n}.fasta \
        -i ${db}.fasta

    sed -i "s/^\(>[^ ]\+\)/\1_c${n}/" ${db}_AA_c${n}.fasta
done


module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/predict_proteomes/select_code.R --prefix ${db}_AA


rm ${db}_AA_c*

chmod 440 ${db}_AA*



### add selected genomes to DB0
cat ${db}.fasta >> DB0.fasta

rm ${db}.fasta




# ------------------------------ other databases ------------------------------ #
for db in 'MGV' 'GPD' 'GVD' 'DEVoC' 'HuVirDB' 'BanfieldLab'
do

    if [ ${db} == 'MGV' ]; then db_fasta=/data/umcg-tifn/DATABASES/MGV/mgv_contigs.fna; fi

    if [ ${db} == 'GPD' ]; then db_fasta=/data/umcg-tifn/DATABASES/GPD/GPD_sequences.fa; fi

    if [ ${db} == 'GVD' ]; then db_fasta=/data/umcg-tifn/DATABASES/GVD/GVDv1_viralpopulations.fna; fi

    if [ ${db} == 'DEVoC' ]; then db_fasta=/data/umcg-tifn/DATABASES/DEVoC/DEVoC.fasta; fi

    if [ ${db} == 'HuVirDB' ]; then db_fasta=/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/HuVirDB-1.0_SafeNames.fasta; fi

    if [ ${db} == 'BanfieldLab' ]; then db_fasta=/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/BanfieldLab.fasta; fi


    sed '1d' /data/umcg-tifn/NL_vir_analysis/get_databases_caudo/${db}_5kb_cMAG_TerL_genomes.txt | cut -d$'\t' -f1 > ${db}.ids

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${db_fasta} ${db}.ids >> DB0.fasta

    rm ${db}.ids

done




# ------------------------------ NL_vir ------------------------------ #
cat /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_genome_fragments.no_rrna.fasta >> DB0.fasta

chmod 440 DB0.fasta
