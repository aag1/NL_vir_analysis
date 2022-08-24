#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=40gb
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-6%2
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ database under consideration ------------------------------ #
cd /data/umcg-tifn/NL_vir_analysis/get_databases_caudo


if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]
then
    db_name='MGV'
    db_fasta=/data/umcg-tifn/DATABASES/MGV/mgv_contigs.fna
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 2 ]
then
    db_name='GPD'
    db_fasta=/data/umcg-tifn/DATABASES/GPD/GPD_sequences.fa
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 3 ]
then
    db_name='GVD'
    db_fasta=/data/umcg-tifn/DATABASES/GVD/GVDv1_viralpopulations.fna
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 4 ]
then
    db_name='DEVoC'
    db_fasta=/data/umcg-tifn/DATABASES/DEVoC/DEVoC.fasta
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 5 ]
then
    db_name='HuVirDB'
    db_fasta=HuVirDB-1.0_SafeNames.fasta
    sed 's/|/_/' /data/umcg-tifn/DATABASES/HuVirDB/HuVirDB-1.0.fasta > ${db_fasta}
    chmod 440 ${db_fasta}
fi

if [ ${SLURM_ARRAY_TASK_ID} -eq 6 ]
then
    db_name='BanfieldLab'
    db_fasta=${db_name}.fasta
    cat /data/umcg-tifn/DATABASES/data_Devoto_2019/15_Lak_phage_genomes.fasta > ${db_fasta}
    cat /data/umcg-tifn/DATABASES/data_Al-Shayeb_2020/349_huge_phage_genomes.fasta >> ${db_fasta}
    cat /data/umcg-tifn/DATABASES/data_Borges_2021_v1/borges_et_al_2021_phage_genomes.fasta >> ${db_fasta}
    chmod 440 ${db_fasta}
fi


echo '-------------------- WORKING WITH '${db_name}' DATABASE --------------------'




# ------------------------------ genomes > 5 kb ------------------------------ #
module purge; module load BBMap; module list

reformat.sh \
    in=${db_fasta} \
    out=${db_name}_5kb.fasta \
    fastawrap=80 \
    minlength=5000




# ------------------------------ genomes with terminal repeats ------------------------------ #
module purge; module load Perl/5.26.1-foss-2018a BioPerl/1.6.924-foss-2018a-Perl-5.26.1; module list

perl /home/umcg-agulyaeva/NL_vir_analysis/taxo_NL_vir/identify_term_rep.pl \
    ${db_name}_5kb.fasta \
    ${db_name}_5kb_cMAG.txt \
    20


sed '1d' ${db_name}_5kb_cMAG.txt | cut -d$'\t' -f1 | sort | uniq > ${db_name}_5kb_cMAG.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    ${db_name}_5kb.fasta \
    ${db_name}_5kb_cMAG.ids > ${db_name}_5kb_cMAG.fasta
    



# ------------------------------ predict proteomes: separate by length ------------------------------ #
module purge; module load BBMap; module list

reformat.sh \
    in=${db_name}_5kb_cMAG.fasta \
    out=${db_name}_5kb_cMAG_under20kb.fasta \
    fastawrap=80 \
    maxlength=19999

reformat.sh \
    in=${db_name}_5kb_cMAG.fasta \
    out=${db_name}_5kb_cMAG_from20kb.fasta \
    fastawrap=80 \
    minlength=20000




# ------------------------------ predict proteomes: genome sequences < 20 kb ------------------------------ #
module purge; module load prodigal; module list

prodigal \
    -p meta \
    -g 11 \
    -f gff \
    -o ${db_name}_5kb_cMAG_under20kb_c11.gff \
    -a ${db_name}_5kb_cMAG_under20kb_c11.fasta \
    -i ${db_name}_5kb_cMAG_under20kb.fasta




# ------------------------------ predict proteomes: genome sequences >= 20 kb ------------------------------ #
i=1
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${db_name}_5kb_cMAG_from20kb_${i}.fasta
        echo $line > $outfile
        i=$((i+1))
    else
        echo $line >> $outfile
    fi
done < ${db_name}_5kb_cMAG_from20kb.fasta


module purge; module load prodigal; module list

for f in ${db_name}_5kb_cMAG_from20kb_*.fasta
do
    b=$(basename ${f} .fasta)

    for n in 11 4 15
    do
        prodigal \
            -p single \
            -g ${n} \
            -f gff \
            -o ${b}_c${n}.gff \
            -a ${b}_c${n}.fasta \
            -i ${f} &
    done

    wait
done




# ------------------------------ predict proteomes: collect data ------------------------------ #
for n in 11 4 15
do
    cat ${db_name}_5kb_cMAG_from20kb_*_c${n}.gff > ${db_name}_5kb_cMAG_c${n}.gff
    cat ${db_name}_5kb_cMAG_from20kb_*_c${n}.fasta > ${db_name}_5kb_cMAG_c${n}.fasta
done


cat ${db_name}_5kb_cMAG_under20kb_c11.gff >> ${db_name}_5kb_cMAG_c11.gff
cat ${db_name}_5kb_cMAG_under20kb_c11.fasta >> ${db_name}_5kb_cMAG_c11.fasta


for n in 11 4 15
do
    sed -i "s/^\(>[^ ]\+\)/\1_c${n}/" ${db_name}_5kb_cMAG_c${n}.fasta
done




# ------------------------------ predict proteomes: select code ------------------------------ #
module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/predict_proteomes/select_code.R \
    --prefix ${db_name}_5kb_cMAG




# ------------------------------ search TerL ------------------------------ #
module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --noali \
    --domtblout markers_vs_${db_name}_5kb_cMAG.txt \
    -o ${db_name}_5kb_cMAG_hmmer_out.txt \
    --cpu ${SLURM_CPUS_PER_TASK} \
    /data/umcg-tifn/NL_vir_analysis/benchmark_markers/markers.hmm \
    ${db_name}_5kb_cMAG_proteomes.fasta


module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/benchmark_markers/select_TerL_hits.R \
    --hmmerF markers_vs_${db_name}_5kb_cMAG.txt \
    --prefix ${db_name}_5kb_cMAG




# ------------------------------ clean up & permissions ------------------------------ #
rm ${db_name}_5kb.fasta
rm ${db_name}_5kb_cMAG.ids
rm ${db_name}_5kb_cMAG.fasta
rm ${db_name}_5kb_cMAG_under20kb*
rm ${db_name}_5kb_cMAG_from20kb*
rm ${db_name}_5kb_cMAG_c*
rm ${db_name}_5kb_cMAG_hmmer_out.txt

chmod 440 markers_vs_${db_name}_5kb_cMAG.txt
chmod 440 ${db_name}_5kb_cMAG*
