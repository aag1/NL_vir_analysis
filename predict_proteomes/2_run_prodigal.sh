#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100%2
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/predict_proteomes

x=$((${SLURM_ARRAY_TASK_ID}-1))

base='NL_vir_genome_fragments_'${x}




# ------------------------------ separate by length ------------------------------ #
module purge; module load BBMap; module list

reformat.sh \
    in=${base}_NT.fasta \
    out=${base}_under20kb.fasta \
    fastawrap=80 \
    maxlength=19999

reformat.sh \
    in=${base}_NT.fasta \
    out=${base}_from20kb.fasta \
    fastawrap=80 \
    minlength=20000




# ------------------------------ genome sequences < 20 kb ------------------------------ #
module purge; module load prodigal; module list

prodigal \
    -p meta \
    -g 11 \
    -f gff \
    -o ${base}_under20kb_c11.gff \
    -a ${base}_under20kb_c11.fasta \
    -i ${base}_under20kb.fasta




# ------------------------------ genome sequences >= 20 kb ------------------------------ #
i=1
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${base}_from20kb_${i}.fasta
        echo $line > $outfile
        i=$((i+1))
    else
        echo $line >> $outfile
    fi
done < ${base}_from20kb.fasta



module purge; module load prodigal; module list

for f in ${base}_from20kb_*.fasta
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
            -i ${f}
    done
done




# ------------------------------ collect data ------------------------------ #
for n in 11 4 15
do
    cat ${base}_from20kb_*_c${n}.gff > ${base}_c${n}.gff
    cat ${base}_from20kb_*_c${n}.fasta > ${base}_c${n}.fasta
done



cat ${base}_under20kb_c11.gff >> ${base}_c11.gff
cat ${base}_under20kb_c11.fasta >> ${base}_c11.fasta



for n in 11 4 15
do
    sed -i "s/^\(>[^ ]\+\)/\1_c${n}/" ${base}_c${n}.fasta
done




# ------------------------------ select code ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/select_code.R --prefix ${base}




# ------------------------------ clean up ------------------------------ #
rm ${base}_under20kb*
rm ${base}_from20kb*
