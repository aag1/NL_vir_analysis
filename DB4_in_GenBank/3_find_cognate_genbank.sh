#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB4_in_GenBank




# ------------------------------ find cognate GenBank genomes ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/find_cognate_genbank.R




# ------------------------------ get publications about the cognate genomes ------------------------------ #
f='DB4_cognate_GenBank_pub_info.txt'

if [ -f ${f} ]; then rm -rf ${f}; fi

echo -e 'gb_id\tauthors\tpub_title\n' > ${f}


ARR=( $(sed '1d' DB4_cognate_GenBank_genera_level.txt | cut -d$'\t' -f4 | sort | uniq) )

for x in "${ARR[@]}"
do

    xml=$(efetch -db nuccore -id ${x} -format xml)

    authors=$(echo ${xml} | xtract -pattern Cit-sub_authors -element Name-std_last | sed 's/\t/, /g')

    pub_title=$(echo ${xml} | xtract -pattern Cit-art_title -element Title_E_name)

    echo -e ${x}'\t'${authors}'\t'${pub_title} >> ${f}

done




# ------------------------------ add publications to the output tables ------------------------------ #
Rscript ${dir}/add_pub_info.R

chmod 440 DB4_cognate_GenBank_*
