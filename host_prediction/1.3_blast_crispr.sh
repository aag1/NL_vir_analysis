#!/bin/bash
#SBATCH --job-name=job1.3
#SBATCH --output=job1.3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/host_prediction




# ------------------------------ phage-host links ------------------------------ #
awk '$13/$14 >= 0.8' blastn_spacers_Shmakov2017_vs_DB3.out > blastn_spacers_Shmakov2017_vs_DB3.spacer80match.out

awk '$13/$14 >= 0.8' blastn_spacers_CRISPRCasdb_vs_DB3.out > blastn_spacers_CRISPRCasdb_vs_DB3.spacer80match.out

module purge; module load R; module list

Rscript ${dir}/crispr_links.R




# ------------------------------ retrieve taxonomy ------------------------------ #
if [ -f DB3_crispr_taxo.txt ]; then rm -rf DB3_crispr_taxo.txt; fi

ARR=( $(sed '1d' DB3_crispr_links.txt | cut -d$'\t' -f2 | sort | uniq) )

for x in "${ARR[@]}"
do
    xml=$(efetch -db nuccore -id ${x} -format xml)

    lineage=$(echo ${xml} | xtract -pattern OrgName -element OrgName_lineage | sed 's/;/,/g')

    organism=$(echo ${xml} | xtract -pattern Org-ref -element Org-ref_taxname)

    genus=$(echo ${xml} | xtract -pattern BinomialOrgName -element BinomialOrgName_genus)

    echo -e ${x}'\t'${lineage}', '${organism}'\t'${genus} >> DB3_crispr_taxo.txt
done




# ------------------------------ link taxonomy ------------------------------ #
Rscript ${dir}/crispr_hosts.R

rm blastn_spacers_*_vs_DB3.spacer80match.out

chmod 440 DB3_crispr_*
