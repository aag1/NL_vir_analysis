#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)

cd /data/umcg-tifn/NL_vir_analysis/test_cenote



module purge; module load Perl/5.26.1-foss-2018a BioPerl/1.6.924-foss-2018a-Perl-5.26.1; module list

perl ${dir}/get_viral_refseq_taxo.pl > viral_refseq_taxo.txt

chmod 440 viral_refseq_taxo.txt



perl ${dir}/get_viral_refseq_aa_nt_ids.pl > viral_refseq_aa_nt_ids.txt

chmod 440 viral_refseq_aa_nt_ids.txt
