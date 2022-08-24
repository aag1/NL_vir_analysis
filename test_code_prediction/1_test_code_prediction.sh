#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/test_code_prediction




# ------------------------------ separate by length ------------------------------ #
module purge; module load BBMap; module list

reformat.sh \
    in=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.fasta \
    out=CRASS_DB_cl_under20kb.fasta \
    fastawrap=80 \
    maxlength=19999

reformat.sh \
    in=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.fasta \
    out=CRASS_DB_cl_from20kb.fasta \
    fastawrap=80 \
    minlength=20000




# ------------------------------ genome sequences < 20 kb ------------------------------ #
module purge; module load prodigal; module list

prodigal \
    -p meta \
    -g 11 \
    -f gff \
    -o CRASS_DB_cl_under20kb_c11.gff \
    -i CRASS_DB_cl_under20kb.fasta




# ------------------------------ genome sequences >= 20 kb ------------------------------ #
# https://www.biostars.org/p/105388/#105464
i=1
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=CRASS_DB_cl_from20kb_${i}.fasta
        echo $line > $outfile
        i=$((i+1))
    else
        echo $line >> $outfile
    fi
done < CRASS_DB_cl_from20kb.fasta



module purge; module load prodigal; module list

for f in CRASS_DB_cl_from20kb_*.fasta
do
    b=$(basename ${f} .fasta)

    for n in 11 4 15
    do
        prodigal \
            -p single \
            -g ${n} \
            -f gff \
            -o ${b}_c${n}.gff \
            -i ${f} &
    done

    wait
done




# ------------------------------ collect data ------------------------------ #
for n in 11 4 15
do
    cat CRASS_DB_cl_from20kb_*_c${n}.gff > CRASS_DB_cl_c${n}.gff
done

cat CRASS_DB_cl_under20kb_c11.gff >> CRASS_DB_cl_c11.gff




# ------------------------------ test code prediction ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/test_code_prediction.R




# ------------------------------ clean up & permissions ------------------------------ #
rm CRASS_DB_cl_from20kb_*.fasta
rm CRASS_DB_cl_from20kb_*_c*.gff
chmod 440 *
