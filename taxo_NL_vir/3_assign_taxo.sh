#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/taxo_NL_vir



perl ${dir}/identify_term_rep.pl \
    /data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.fasta \
    NL_vir_genome_fragments_circ.txt \
    20



module purge; module load R; module list

Rscript /home/umcg-agulyaeva/NL_vir_analysis/benchmark_markers/select_TerL_hits.R \
    --hmmerF markers_vs_NL_vir.txt \
    --prefix NL_vir

Rscript ${dir}/assign_taxo.R



chmod 440 NL_vir_genome_fragments_circ.txt
chmod 440 NL_vir_TerL_*.txt
chmod 440 NL_vir_genome_fragments.no_rrna.taxo.txt
