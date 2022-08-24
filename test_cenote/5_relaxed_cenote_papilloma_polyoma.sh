#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/test_cenote




# ------------------------------ extract RefSeq Papillomaviridae & Polyomaviridae ------------------------------ #
grep 'Papillomaviridae' viral_refseq_taxo.txt | cut -d$'\t' -f1 > papilloma_polyoma.ids
grep 'Polyomaviridae' viral_refseq_taxo.txt | cut -d$'\t' -f1 >> papilloma_polyoma.ids
grep 'Hepadnaviridae' viral_refseq_taxo.txt | cut -d$'\t' -f1 >> papilloma_polyoma.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    /data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_genomic.fna \
    papilloma_polyoma.ids > papilloma_polyoma.fasta




# ------------------------------ RefSeq Papillomaviridae & Polyomaviridae: circularized? ------------------------------ #
perl /home/umcg-agulyaeva/NL_vir_analysis/taxo_NL_vir/identify_term_rep.pl \
    papilloma_polyoma.fasta \
    papilloma_polyoma_circular.txt \
    4




# ------------------------------ run relaxed cenote ------------------------------ #
module purge; module load Anaconda3; module list
source activate cenote-taker2_env; conda list

python /data/umcg-agulyaeva/SOFTWARE/Cenote-Taker2/unlimited_breadsticks.py \
    --contigs papilloma_polyoma.fasta \
    --run_title papilloma_polyoma \
    --virus_domain_db 'virion' \
    --minimum_length_circular 3000 \
    --minimum_length_linear 3000 \
    --circ_minimum_hallmark_genes 1 \
    --lin_minimum_hallmark_genes 1 \
    --prune_prophage True \
    --filter_out_plasmids True \
    --mem 4 \
    --cpu ${SLURM_CPUS_PER_TASK}

conda deactivate

mv papilloma_polyoma/papilloma_polyoma_CONTIG_SUMMARY.tsv .
mv papilloma_polyoma/papilloma_polyoma_PRUNING_INFO_TABLE.tsv .

chmod 440 papilloma_polyoma_circular.txt
chmod 440 papilloma_polyoma_CONTIG_SUMMARY.tsv
chmod 440 papilloma_polyoma_PRUNING_INFO_TABLE.tsv

rm -rf papilloma_polyoma
rm papilloma_polyoma.ids
rm papilloma_polyoma.fasta
rm papilloma_polyoma.over_3000nt.fasta




# ------------------------------ summarize results of the relaxed cenote ------------------------------ #
module purge; module load R; module list

Rscript ${dir}/summarize_results_relaxed_cenote.R
