#!/bin/bash
#SBATCH --job-name=job2.3
#SBATCH --output=job2.3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/host_prediction




# ------------------------------ select hits ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/select_host_fragm_hits.R




# ------------------------------ retrieve taxonomy ------------------------------ #
if [ -f DB3_host_fragm_taxo.txt ]; then rm -rf DB3_host_fragm_taxo.txt; fi

ARR=( $(sed '1d' DB3_host_fragm_sele_hits.txt | cut -d$'\t' -f3 | sort | uniq) )

for x in "${ARR[@]}"
do
    xml=$(efetch -db nuccore -id ${x} -format xml)

    lineage=$(echo ${xml} | xtract -pattern OrgName -element OrgName_lineage | sed 's/;/,/g')

    organism=$(echo ${xml} | xtract -pattern Org-ref -element Org-ref_taxname)

    genus=$(echo ${xml} | xtract -pattern BinomialOrgName -element BinomialOrgName_genus)

    echo -e ${x}'\t'${lineage}', '${organism}'\t'${genus} >> DB3_host_fragm_taxo.txt
done




# ------------------------------ link taxonomy ------------------------------ #
Rscript ${dir}/assign_host_fragm_hosts.R

chmod 440 DB3_host_fragm_*
