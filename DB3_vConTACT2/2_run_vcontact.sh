#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=120gb
#SBATCH --partition=himem
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/DB3_vConTACT2


echo 'protein_id,contig_id,keywords' > gene_to_genome.csv
awk 'BEGIN{FS="\t"; OFS=","} (NR>1) {print $3, $2, "None_provided"}' proteomes.txt >> gene_to_genome.csv


module purge; module load Anaconda3; module list
source activate vContact2; conda list

vcontact2 \
	--raw-proteins proteomes.fasta \
	--proteins-fp gene_to_genome.csv \
	--db 'None' \
	--output-dir 'vcontact_output' \
	--c1-bin '/home/umcg-agulyaeva/SOFTWARE/cluster_one-1.0.jar' \
	--threads ${SLURM_CPUS_PER_TASK}

conda deactivate


mv vcontact_output/c1.ntw .
mv vcontact_output/viral_cluster_overview.csv .
mv vcontact_output/genome_by_genome_overview.csv .
rm -rf vcontact_output
chmod 440 gene_to_genome.csv c1.ntw viral_cluster_overview.csv genome_by_genome_overview.csv
