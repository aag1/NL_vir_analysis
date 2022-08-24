#!/bin/bash
#SBATCH --job-name=job2.1
#SBATCH --output=job2.1_%A.out
#SBATCH --mem=40gb
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/NL_vir_analysis/cluster_genomes
base='DB0'
### see https://bitbucket.org/berkeleylab/checkv/src/master/



### create a blast+ database
module purge; module load BLAST+; module list

makeblastdb \
    -in ${base}.fasta \
    -dbtype nucl \
    -out ${base}



### perform all-vs-all blastn of sequences
blastn \
    -query ${base}.fasta \
    -db ${base} \
    -max_target_seqs 10000 \
    -out ${base}_blast.tsv \
    -outfmt '6 std qlen slen' \
    -num_threads ${SLURM_CPUS_PER_TASK}



rm ${base}.n*
chmod 440 ${base}_blast.tsv
