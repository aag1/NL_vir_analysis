#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/DB3_TerL_tree




module purge; module load R; module list

Rscript ${dir}/summary_data4plot.R

Rscript ${dir}/summary_phage_families.R

Rscript ${dir}/summary_host_phyla.R

Rscript ${dir}/summary_pheno_assoc.R

Rscript ${dir}/plot_trees.R




chmod 440 DB3_data4plot.txt
chmod 440 DB3_phage_families_matrix.txt DB3_phage_families_extended_via_mrca.txt
chmod 440 DB3_predicted_host_phyla.txt
chmod 440 DB4_pheno_assoc_pch.txt
chmod 440 DB*_TerL_tree.pdf
