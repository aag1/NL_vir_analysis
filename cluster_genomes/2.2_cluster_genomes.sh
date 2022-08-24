#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A.out
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/cluster_genomes
base='DB0'
### see https://bitbucket.org/berkeleylab/checkv/src/master/



### calculate pairwise ANI by combining local alignments between sequence pairs
module purge; module load Anaconda3; module list

source activate CheckV
conda list

python /home/umcg-agulyaeva/SOFTWARE/CheckV_scripts/anicalc.py \
    -i ${base}_blast.tsv \
    -o ${base}_ani.tsv



### perform UCLUST-like clustering
python /home/umcg-agulyaeva/SOFTWARE/CheckV_scripts/aniclust.py \
    --fna ${base}.fasta \
    --ani ${base}_ani.tsv \
    --out ${base}_clusters_ANI95_AF85.tsv \
    --min_ani 95 \
    --min_qcov 0 \
    --min_tcov 85

conda deactivate



chmod 440 ${base}_ani.tsv ${base}_clusters_*.tsv
