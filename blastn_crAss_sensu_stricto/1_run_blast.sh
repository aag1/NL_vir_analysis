#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



for id in 'NC_024711.1' 'MGV-GENOME-0359371'
do
    echo ${id} > ${id}.id

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB0.fasta \
        ${id}.id > ${id}.fa
done



/home/umcg-agulyaeva/SOFTWARE/ncbi-blast-2.12.0+/bin/blastn \
    -task 'blastn' \
    -query NC_024711.1.fa \
    -subject MGV-GENOME-0359371.fa \
    -out crAss_sensu_stricto_vs_MGV-GENOME-0359371.txt \
    -evalue 0.05 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore nident'



module purge; module load R; module list

Rscript qcov_pct.R



rm *id *fa
