sessionInfo()




# ------------------------------ NL_vir contigs info ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/taxo_NL_vir/NL_vir_genome_fragments.no_rrna.taxo.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)



### total number & % from different cohorts
nrow(t)
# 58776

for (k in unique(t$cohort)) { cat(k, round(sum(t$cohort == k) / nrow(t) * 100), '%\n') }
# LLD 45 %
# LLD2 21 %
# 300OB 15 %
# IBD 19 %



### NL_vir contigs with terminal repeats (NB! 'circularized' is not a suitable term to designate a genome with terminal repeats!)
table(t$contig_fragm_circ)
#            DTR DTR;ITR     ITR
#  57163    1567       1      45

sum(t$contig_fragm_circ != '')
# 1613

round(sum(grepl('DTR', t$contig_fragm_circ)) / sum(t$contig_fragm_circ != '') * 100, 0)
# 97

round(sum(grepl('ITR', t$contig_fragm_circ)) / sum(t$contig_fragm_circ != '') * 100, 0)
# 3



### NL_vir contigs taxonomy
table(t$contig_taxo)
#                  Adenoviridae Caudoviricetes
#          19023              1          39752

t[t$contig_taxo == 'Adenoviridae', c('final_name', 'contig_fragm_len', 'contig_fragm_circ')]
#   final_name contig_fragm_len contig_fragm_circ
# NL_vir016001            17402



### prophage NL_vir contigs
p <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

sum(p$prophage_cenote[ p$source == 'NL_vir' ])
# 15570



### NL_vir contigs with alternative genetic codes
table(p$genome_code[ p$source == 'NL_vir' ])
#   c11   c15    c4
# 53049  5706    21



### rRNA NL_vir contigs
# wc -l /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_rrna_contigs.ids
# 8




# ------------------------------ number of sequences from databases ------------------------------ #
# cd /data/umcg-tifn/NL_vir_analysis/get_databases_caudo/; wc -l *_5kb_cMAG_TerL_genomes.txt
#    459 BanfieldLab_5kb_cMAG_TerL_genomes.txt
#    138 DEVoC_5kb_cMAG_TerL_genomes.txt
#  10871 GPD_5kb_cMAG_TerL_genomes.txt
#    937 GVD_5kb_cMAG_TerL_genomes.txt
#   1661 HuVirDB_5kb_cMAG_TerL_genomes.txt
#  19695 MGV_5kb_cMAG_TerL_genomes.txt

# cd /data/umcg-tifn/NL_vir_analysis/cluster_genomes/; wc -l RefSeq.ids; wc -l Benler_2020.ids
# 6049 RefSeq.ids
# 1480 Benler_2020.ids




# ------------------------------ number of DB1 vOTUs ------------------------------ #
tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, stringsAsFactors = F)

nrow(tab)
# 100060

length(unique(tab$sp_repres))
# 30461




# ------------------------------ DB2 info ------------------------------ #
DF <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)

TAB <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_circ_PlasX_score.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)



### number of vOTUs
nrow(DF)
# 15196



### number of caudo vOTUs represented by genomes with terminal repeats
sum((DF$simple_termini == 'TR') & (DF$simple_taxo == 'Caudoviricetes'))
# 2106



### taxonomic composition
table(DF$simple_taxo)
#    Adenoviridae   Caudoviricetes     Microviridae Papillomaviridae
#               2            10514                2                2
#  Polyomaviridae       unassigned
#               1             4675

table(DF$simple_taxo) / nrow(DF) * 100
#    Adenoviridae   Caudoviricetes     Microviridae Papillomaviridae
#     0.013161358     69.189260332      0.013161358      0.013161358
#  Polyomaviridae       unassigned
#     0.006580679     30.764674914



### Are plasmids more prevalent among DTR genomes representing vOTUs that were classified as Caudoviricetes, or that remained taxonomically unassigned?
DF$DTR <- sapply(rownames(DF), function (x) grepl('DTR', tab$genome_circ[tab$genome_id == x]))

df <- DF[DF$DTR, ]

df$plasmid <- ifelse(rownames(df) %in% TAB$contig[ TAB$score > 0.9 ], 1, 0)

for (x in c('Caudoviricetes', 'Unassigned')) {

    n <- sum((df$simple_taxo == x) & (df$plasmid == 1)) / sum(df$simple_taxo == x) * 100
    n <- round(n, 2)
    cat('% plasmids among DTR genomes representing vOTUs classified as', x, ':', n, '%\n')

}

# % plasmids among DTR genomes representing vOTUs classified as Caudoviricetes : 2.04 %
# % plasmids among DTR genomes representing vOTUs classified as Unassigned : 20 %




# ------------------------------ Top 10 DB3 phages most prevalent in DEVoC HA ------------------------------ #
tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, stringsAsFactors = F)

d1 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_TerL_tree/DB3_data4plot.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

d2 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_TerL_tree/DB3_phage_families_extended_via_mrca.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]



d1 <- d1[, grep('^p_samples_', colnames(d1), value = T)]

d1 <- d1[order(d1$p_samples_DEVoC_HA, decreasing = T), ]

d1[, colnames(d2)] <- d2[rownames(d1), ]

d1$db4 <- ifelse(rownames(d1) %in% ids4, 1, 0)

print(d1[1:10, ], width = 500)

#                                               p_samples_LLD p_samples_LLD2 p_samples_300OB p_samples_IBD p_samples_DEVoC_all p_samples_DEVoC_HA    phage_family phage_family_ext db4
# uvig_191910                                      12.6872247      20.118343       13.758389     6.1538462           15.748031          28.846154      crAss-like       crAss-like   1
# MGV-GENOME-0193745                                0.5286344       0.000000        0.000000     0.0000000           20.866142          21.153846                                    0
# MGV-GENOME-0279285                                8.7224670      16.272189       23.489933     6.5934066            5.511811          19.230769                                    1
# SRS415869_NODE_11_length_123658_cov_86.751276     5.4625551      12.721893        7.718121     1.9780220            5.905512          17.307692                                    1
# uvig_249048                                       1.5859031       2.662722        1.006711     0.4395604            9.055118          15.384615      crAss-like       crAss-like   0
# uvig_205988                                       6.4317181      11.834320       10.402685     2.6373626            8.267717          15.384615                                    1
# uvig_535382                                       1.1453744       2.366864        2.348993     0.2197802            2.755906          13.461538                                    0
# MGV-GENOME-0289925                                4.4052863       7.692308        9.060403     1.5384615            2.362205           9.615385                                    1
# MGV-GENOME-0359951                                1.8502203       3.254438        2.013423     0.8791209            5.118110           9.615385      crAss-like       crAss-like   0
# MGV-GENOME-0351725                                2.8193833       7.692308        1.342282     1.9780220            5.118110           9.615385 Flandersviridae  Flandersviridae   1



print(tab[tab$sp_repres == 'MGV-GENOME-0193745', c('genome_id', 'taxo_refseq')], width = 500)

#                genome_id                                                                                                                                       taxo_refseq
# 4087         NC_024390.1    Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Siphoviridae;Mccleskeyvirinae;Limdunavirus;Leuconostoc phage LN03
# 4088         NC_020870.1    Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Siphoviridae;Mccleskeyvirinae;Limdunavirus;Leuconostoc phage LN04
# 4089         NC_024385.1 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Siphoviridae;Mccleskeyvirinae;Limdunavirus;Leuconostoc phage phiLN12
# 18190       NL_vir012145
# 66570 MGV-GENOME-0193745




# ------------------------------ DB3 composition numbers ------------------------------ #
d <- read.table('/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

ids3 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]



### DB3 constitutes ? % of DB2
round(length(ids3) / nrow(d) * 100)
# 12



### DB3 constitutes ? % of NL_vir virus-like contigs
# grep -c '>' /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/NL_vir_from_DB3_vOTUs.fasta
# 17193

# grep -c '>' /data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_genome_fragments.no_rrna.fasta
# 58776

round(17193 / 58776 * 100)
# 29



### DB3 vOTU sources
round(table(d[ids3, 'simple_source']) / length(ids3) * 100)
#      Both  Databases This study
#        56         25         19




# ------------------------------ DB3 detection numbers ------------------------------ #
for (x in c('LLD', 'LLD2', '300OB', 'IBD', 'DEVoC_HA', 'DEVoC_all')) {

    if (x %in% c('DEVoC_HA', 'DEVoC_all')) {

        f <- '/data/umcg-tifn/NL_vir_analysis/DEVoC_reads/DEVoC_DB3_mapping_summary.txt'

    } else {

        f <- paste0('/data/umcg-tifn/NL_vir_analysis/map_reads/', x, '_DB3_mapping_summary.txt')

    }


    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)


    if (x == 'DEVoC_HA') {

       d <- read.table('/home/umcg-agulyaeva/NL_vir_analysis/DEVoC_reads/DEVoC_adult_cohort_PRJNA722819.tsv', sep = '\t', header = T, stringsAsFactors = F)
       sele <- d$run_accession[ grepl('^HA', d$sample_alias) ]
       t <- t[t$sample_id %in% sele, ]

    }


    cat(x, 'cohort:\n')
    cat('# vOTUs detected mean', round(mean(t$num_sp_cov75)), '\n')
    cat('% reads mapped mean', round(mean(t$pct_reads_mapped), 2), '%\n\n')

}

# LLD cohort:
# # vOTUs detected mean 7 
# % reads mapped mean 0.68 %

# LLD2 cohort:
# # vOTUs detected mean 10 
# % reads mapped mean 0.75 %

# 300OB cohort:
# # vOTUs detected mean 9 
# % reads mapped mean 0.74 %

# IBD cohort:
# # vOTUs detected mean 5 
# % reads mapped mean 0.6 %

# DEVoC_HA cohort:
# # vOTUs detected mean 6 
# % reads mapped mean 6.92 %

# DEVoC_all cohort:
# vOTUs detected mean 5 
# % reads mapped mean 7.88 %




# ------------------------------ DB4 genomes with taxonomic info ------------------------------ #
tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, stringsAsFactors = F)

ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]



### RefSeq genomes in vOTUs
for (x in ids4) {

    idx <- which((tab$sp_repres == x) & (tab$source == 'RefSeq'))

    if (length(idx) > 0) { print(tab[idx, c('sp_repres', 'genome_id', 'taxo_refseq')], width = 500) }

}

#               sp_repres   genome_id                                                                                                                                                           taxo_refseq
# 5237 MGV-GENOME-0305083 NC_047911.1 Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Myoviridae;Lagaffevirus;Faecalibacterium virus Lagaffe;Faecalibacterium phage FP_Lagaffe



### BanfieldLab genomes in vOTUs
for (x in ids4) {

    idx <- which((tab$sp_repres == x) & (tab$source == 'BanfieldLab'))

    if (length(idx) > 0) { print(tab[idx, c('sp_repres', 'genome_id', 'taxo_BanfieldLab')], width = 500) }

}

# nothing



### Benler_2020 genomes in vOTUs
for (x in ids4) {

    idx <- which((tab$sp_repres == x) & (tab$source == 'Benler_2020') & (tab$taxo_Benler_2020 != 'Uroviricota,Caudoviricetes,Caudovirales,,,'))

    if (length(idx) > 0) { print(tab[idx, c('sp_repres', 'genome_id', 'taxo_Benler_2020')], width = 500) }

}

#                sp_repres      genome_id                                          taxo_Benler_2020
# 64881 MGV-GENOME-0351725 OJML01000036.1 Uroviricota,Caudoviricetes,Caudovirales,Flandersviridae,,

#         sp_repres      genome_id                           taxo_Benler_2020
# 66179 uvig_191910 OLGH01000138.1 Uroviricota,Caudoviricetes,Crassvirales,,,

#                sp_repres      genome_id                                         taxo_Benler_2020
# 64909 MGV-GENOME-0328053 OLVO01000312.1 Uroviricota,Caudoviricetes,Caudovirales,group4986,,I4986




# ------------------------------ DB4 crAss-like vOTUs ------------------------------ #
tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, stringsAsFactors = F)

D1 <- read.table('/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt', sep = '\t', row.names = 2, header = T, stringsAsFactors = F)

D2 <- read.table('/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)



### MGV-GENOME-0359371
t <- tab[tab$sp_repres == 'MGV-GENOME-0359371', ]


table(t$source)
#    GPD HuVirDB     MGV  NL_vir
#     20       9      43      20


v <- t$taxo_our_crAss
v <- v[v != '']
v <- unique(D1[v, 'cl_repres'])
unique(D2[v, 'genus'])
# alpha17



### uvig_191910
t <- tab[tab$sp_repres == 'uvig_191910', ]


table(t$source)
# Benler_2020         GPD     HuVirDB         MGV      NL_vir
#           1          97           4         207         178


v <- t$taxo_our_crAss
v <- v[v != '']
v <- unique(D1[v, 'cl_repres'])
unique(D2[v, 'genus'])
# gamma4




# ------------------------------ Supressor tRNAs ------------------------------ #
### DB3
ids3 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]

tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

R <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_tRNA_scan/DB3_tRNAscan.txt', sep = '\t', header = T, stringsAsFactors = F)
R$genome_id <- sub(' $', '', R$genome_id)


DF <- data.frame(
    genome_id = ids3,
    genome_code = tab[ids3, 'genome_code'],
    sup_tRNA_anticodon = sapply(ids3, function (x) {

        idx <- which((R$genome_id == x) & (R$tRNA_type == 'Sup'))

        if (length(idx) > 0) {
            paste(unique(R$anticodon[idx]), collapse = ';')
        } else { '' }

    }),
    stringsAsFactors = F
)


X <- aggregate(DF$genome_id, by = DF[, c('genome_code', 'sup_tRNA_anticodon')], FUN = length)
colnames(X) <- c('genome_code', 'sup_tRNA_anticodon', 'db3_genomes_number')
X
#   genome_code sup_tRNA_anticodon db3_genomes_number
#1          c11                                  1581
#2          c15                                   248
#3           c4                                     5
#4          c11                CTA                 13
#5          c15                CTA                 36
#6          c15            CTA;TTA                  1
#7           c4                TCA                  1
#8          c11            TCA;CTA                  1
#9          c11        TCA;CTA;TTA                  1
#10         c11                TTA                 12


sum(X$db3_genomes_number[(X$genome_code == 'c15') & grepl('CTA', X$sup_tRNA_anticodon)]) / sum(X$db3_genomes_number[(X$genome_code == 'c15')]) * 100
# 12.98246



### crAss378
tab <- read.table('/data/umcg-tifn/NL_vir_analysis/test_code_prediction/crAss378_predicted_code.txt', sep = '\t', header = T, stringsAsFactors = F)

R <- read.table('/data/umcg-tifn/NL_vir_analysis/test_code_prediction/crAss378_tRNAscan.txt', sep = '\t', header = T, stringsAsFactors = F)
R$genome_id <- sub(' $', '', R$genome_id)


DF <- tab[, c('genome_id', 'sele_code', 'real_code')]

DF$sup_tRNA_anticodon = sapply(DF$genome_id, function (x) {

    idx <- which((R$genome_id == x) & (R$tRNA_type == 'Sup'))

    if (length(idx) > 0) {
        paste(unique(R$anticodon[idx]), collapse = ';')
    } else { '' }

})


X <- aggregate(DF$genome_id, by = DF[, c('sele_code', 'real_code', 'sup_tRNA_anticodon')], FUN = length)
colnames(X) <- c('sele_code', 'real_code', 'sup_tRNA_anticodon', 'crAss378_genomes_number')
X
#   sele_code real_code sup_tRNA_anticodon crAss378_genomes_number
#1        c11       c11                                        241
#2         c4       c11                                          1
#3        c11       c15                                          5
#4        c15       c15                                         46
#5         c4        c4                                          8
#6        c11       c11                CTA                       2
#7        c11       c15                CTA                       5
#8        c15       c15                CTA                      66
#9         c4        c4                TCA                       2
#10       c11        c4        TCA;CTA;TTA                       1
#11       c11       c11                TTA                       1




# ------------------------------ DB3 vOTUs with RefSeq sequences ------------------------------ #
#awk -F'\t' '(NR > 1) && ($11 != "")' /data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3_table1.txt | cut -d$'\t' -f1 | sort | uniq | wc -l
# 25




# ------------------------------ pheno assoc ------------------------------ #
#sed '1d' /data/umcg-tifn/NL_vir_analysis/host_prediction/DB4_compare_host_predictions.txt | cut -d$'\t' -f1 | sort | uniq | wc -l
# 32

(32 - 3) / 32 * 100
# 90.625



for (x in c('300OB', 'IBD')) {

    cat('\n', x, 'cohort pheno assoc:\n\n')


    t <- read.table(
        paste0('/home/umcg-agulyaeva/NL_vir_analysis/DB3_TerL_tree/from_Gearshift/LLD_vs_', x, '.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )


    for (y in c('unadjusted', 'adjusted')) {

        sele_fdr <- (t[, paste0('FDR.glm', ifelse(y == 'unadjusted', 1, 2))] < 0.05)

        for (z in c('overrepresented', 'underrepresented')) {

            if (z == 'overrepresented')  { sele_preval <- (t[, 'LLD_positive_samples'] < t[, paste0(sub('^3', 'X3', x), '_positive_samples')]) }
            if (z == 'underrepresented') { sele_preval <- (t[, 'LLD_positive_samples'] > t[, paste0(sub('^3', 'X3', x), '_positive_samples')]) }

            N <- sum(sele_fdr & sele_preval)
            cat(z, 'in', x, 'and', y, 'for host:', N, '\n\n')

        }

    }

}

# 300OB cohort pheno assoc:

# overrepresented in 300OB and unadjusted for host: 8 

# underrepresented in 300OB and unadjusted for host: 0 

# overrepresented in 300OB and adjusted for host: 2 

# underrepresented in 300OB and adjusted for host: 0 


# IBD cohort pheno assoc:

# overrepresented in IBD and unadjusted for host: 7 

# underrepresented in IBD and unadjusted for host: 16 

# overrepresented in IBD and adjusted for host: 3 

# underrepresented in IBD and adjusted for host: 5  
