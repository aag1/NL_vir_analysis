#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=1gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




cd /data/umcg-tifn/NL_vir_analysis/DB3_vConTACT2

if [ -f proteomes.fasta ]; then rm -rf proteomes.fasta; fi

if [ -f proteomes.txt ]; then rm -rf proteomes.txt; fi

echo -e "taxonomy\tgenome\tprotein" > proteomes.txt




# ------------------------------ phylum Mirusviricota proteomes ------------------------------ #
file_AA=/data/umcg-tifn/DATABASES/data_Gaia_2022/GOEV_Genes_AA_fa



if [ -f Mirusviricota_111_AA.ids ]; then rm -rf Mirusviricota_111_AA.ids; fi

arr=( $(grep '>' /data/umcg-tifn/DATABASES/data_Gaia_2022/Mirus_111_MAGs.fa | sed 's/^>//') )

for x in "${arr[@]}"
do

    grep "^>[0-9]\+ contig:${x};" ${file_AA} | sed 's/^>//' | cut -d' ' -f1 >> Mirusviricota_111_AA.ids

done



/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${file_AA} Mirusviricota_111_AA.ids > proteomes.fasta

sed -i 's/^>\([0-9]\+\) contig:\([A-Z0-9_]\+\);.\+$/>\2_\1/' proteomes.fasta

rm Mirusviricota_111_AA.ids



grep '>' proteomes.fasta | sed 's/^>\(TARA_[A-Z]\+_NCLDV_[0-9]\+\)\(_[0-9]\+_[0-9]\+\)$/Mirusviricota\t\1\t\1\2/' >> proteomes.txt




# ------------------------------ class Herviviricetes & Caudoviricetes proteomes ------------------------------ #
file1=/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_taxo.txt

file2=/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt

file3=/data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_protein.faa



for x in 'Herviviricetes' 'Caudoviricetes'
do

    if [ -f ${x}_AA.ids ]; then rm -rf ${x}_AA.ids; fi


    genomes=( $(awk -v var=${x} -F'\t' '{if ($7 ~ "\\<"var"\\>" && $2 ~ /complete genome/) print $1}' ${file1}) )

    for g in "${genomes[@]}"
    do

        proteins=( $(awk -v var=${g} -F'\t' '{if ($1 == var) print $2}' ${file2}) )

        for p in "${proteins[@]}"
        do

            echo ${p} >> ${x}_AA.ids

            echo -e "${x}\t${g}\t${p}" >> proteomes.txt

        done

    done


    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${file3} ${x}_AA.ids >> proteomes.fasta

    rm ${x}_AA.ids

done




# ------------------------------ proteomes of DB3 vOTU representatives ------------------------------ #
genomes=( $(cat /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.ids) )

for g in "${genomes[@]}"
do
    s=$(awk -v var=${g} '{if ($1 == var) print $8}' /data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt)


    if [ ${s} == 'RefSeq' ]; then

        if cut -d$'\t' -f2 proteomes.txt | grep -q -x ${g}
        then
            echo 'RefSeq genome '$g' is already in proteomes.txt!'
        else
            echo 'RefSeq genome '$g' is missing from proteomes.txt!'
        fi

        continue

    elif [ ${s} == 'NL_vir' ]; then

        f1='/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt'
        f2='/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.fasta'

    elif [ ${s} == 'Benler_2020' ]; then

        f1='/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.txt'
        f2='/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.fasta'

    else

        f1='/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/'${s}'_5kb_cMAG_proteomes.txt'
        f2='/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/'${s}'_5kb_cMAG_proteomes.fasta'

    fi


    if [ -f ${g}_AA.ids ]; then rm -rf ${g}_AA.ids; fi

    proteins=( $(awk -F'\t' -v var=${g} '{if ($1 == var) print $2}' ${f1}) )

    for p in "${proteins[@]}"
    do

        echo ${p} >> ${g}_AA.ids

        echo -e "DB3\t${g}\t${p}" >> proteomes.txt

    done

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${f2} ${g}_AA.ids >> proteomes.fasta

    rm ${g}_AA.ids

done




chmod 440 proteomes.fasta proteomes.txt
