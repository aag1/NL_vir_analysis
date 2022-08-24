sessionInfo()




# ------------------------------ function ------------------------------ #
term_rep_type <- function (seq_id, seq_len, t) {

    idx <- which(t$seq_id == seq_id)

    if (length(idx) == 0) {

        res <- ''

    } else {

        res <- paste(t$term_rep_type[idx], collapse = ';')

    }


    idx <- which((t$seq_id == seq_id) & (t$term_rep_type == 'ITR'))

    if (length(idx) == 1) {

        if (t$term_rep_len[idx] == seq_len) {

            cat(seq_id, 'is identical to its reverse complement!\n')

        }

    }


    return(res)

}




# ------------------------------ Viral RefSeq ------------------------------ #
tab <- read.table(
    'DB0_lengths_GC.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)



v <- read.table(
    'RefSeq.ids',
    header = F,
    sep = '\t',
    stringsAsFactors = F
)[, 1]



t1 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_taxo.txt',
    header = T,
    sep = '\t',
    quote = '',
    fill = T,
    row.names = 1,
    stringsAsFactors = F
)



t2 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/benchmark_markers/viral_refseq_TerL_genomes.txt',
    header = T,
    sep = '\t',
    row.names = 1,
    stringsAsFactors = F
)



t3 <- read.table(
    'RefSeq_AA_SeleCode.txt',
    header = T,
    sep = '\t',
    row.names = 1,
    stringsAsFactors = F
)



w <- read.table(
    'RefSeq_cMAG.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = F
)



k <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt',
    header = T,
    sep = '\t',
    quote = '',
    stringsAsFactors = F
)



DF <- data.frame(
    genome_id        = v,
    genome_len       = tab[v, 'Length'],
    genome_pgc       = tab[v, 'X.GC'],
    genome_code      = t3[v, 'sele_code'],
    genome_circ      = sapply(v, function (x) { term_rep_type(x, tab[x, 'Length'], w) }),
    prophage_cenote  = NA,
    cohort           = NA,
    source           = 'RefSeq',
    source_formal    = 'Viral RefSeq 209',
    terL_id          = sapply(v, function (x) ifelse(x %in% rownames(t2), t2[x, 'protein_id'], '')),
    terL_coo         = '',
    taxo_refseq      = t1[v, 'taxonomy'],
    taxo_Benler_2020 = '',
    taxo_BanfieldLab = '',
    taxo_our_crAss   = '',
    stringsAsFactors = F
)



sele <- which(DF$terL_id != '')
DF$terL_coo[sele] <- sapply(DF$terL_id[sele], function (x) { k$coded_by[k$protein_id == x] })




# ------------------------------ NL_vir ------------------------------ #
t1 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/taxo_NL_vir/NL_vir_genome_fragments.no_rrna.taxo.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = F
)

t1$contig_len <- as.numeric(sub('^NODE_[0-9]+_length_([0-9]+)_cov_[0-9\\.]+$', '\\1', t1$contig))

t1$prophage_cenote <- ifelse(t1$contig_len == t1$contig_fragm_len, 0, 1)



t2 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_SeleCode.txt',
    header = T,
    sep = '\t',
    row.names = 1,
    stringsAsFactors = F
)



t3 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/taxo_NL_vir/NL_vir_TerL_genomes.txt',
    header = T,
    sep = '\t',
    row.names = 1,
    stringsAsFactors = F
)



k <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = F
)



df <- data.frame(
    genome_id        = t1$final_name,
    genome_len       = t1$contig_fragm_len,
    genome_pgc       = tab[t1$final_name, 'X.GC'],
    genome_code      = t2[t1$final_name, 'sele_code'],
    genome_circ      = t1$contig_fragm_circ,
    prophage_cenote  = t1$prophage_cenote,
    cohort           = t1$cohort,
    source           = 'NL_vir',
    source_formal    = paste0('This study, ', sub('^LLD2$', 'LLD follow-up', t1$cohort), ' cohort'),
    terL_id          = sapply(t1$final_name, function (x) ifelse(x %in% rownames(t3), t3[x, 'protein_id'], '')),
    terL_coo         = '',
    taxo_refseq      = '',
    taxo_Benler_2020  = '',
    taxo_BanfieldLab = '',
    taxo_our_crAss   = t1$crass_paper,
    stringsAsFactors = F
)



sele <- which(df$terL_id != '')
df$terL_coo[sele] <- sapply(df$terL_id[sele], function (x) {
                        i <- which(k$protein_id == x)
                        s <- ifelse(k$strain[i] == 1, 'f', 'r')
                        paste0(s, ';', k$start[i], ';', k$end[i])
})



DF <- rbind(DF, df)




# ------------------------------ Benler et al. 2020 ------------------------------ #
t1 <- read.table(
    '/data/umcg-tifn/DATABASES/data_Benler_2020/40168_2021_1017_MOESM3_ESM.csv',
    header = T,
    sep = ',',
    quote = '',
    fill = T,
    stringsAsFactors = F
)

t1 <- t1[t1$Phylum == 'Uroviricota', ]

t1$Taxonomy <- apply(t1[, c('Phylum', 'Class', 'Order', 'Family', 'Subfamily', 'Genus')], 1, function (v) paste(v, collapse = ','))



t2 <- read.table(
    'Benler_2020_on_Fig1_TerL.txt',
    header = F,
    sep = '\t',
    stringsAsFactors = F
)

v <- setNames(t2[, 2], t2[, 1])



t3 <- read.table(
    'Benler_2020_AA_SeleCode.txt',
    header = T,
    sep = '\t',
    row.names = 1,
    stringsAsFactors = F
)



w <- read.table(
    'Benler_2020_cMAG.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = F
)



k <- read.table(
    '/data/umcg-tifn/DATABASES/data_Benler_2020/gut_phage_annotations.csv',
    header = T,
    sep = ',',
    quote = '',
    fill = T,
    stringsAsFactors = F
)



df <- data.frame(
    genome_id        = t1$representative,
    genome_len       = t1$length,
    genome_pgc       = tab[t1$representative, 'X.GC'],
    genome_code      = t3[t1$representative, 'sele_code'],
    genome_circ      = sapply(1:nrow(t1), function (i) { term_rep_type(t1$representative[i], t1$length[i], w) }),
    prophage_cenote  = NA,
    cohort           = NA,
    source           = 'Benler_2020',
    source_formal    = 'Benler et al. 2021',
    terL_id          = sapply(t1$representative, function (x) ifelse(x %in% names(v), v[x], '')),
    terL_coo         = '',
    taxo_refseq      = '',
    taxo_Benler_2020  = t1$Taxonomy,
    taxo_BanfieldLab = '',
    taxo_our_crAss   = '',
    stringsAsFactors = F
)



sele <- which(df$terL_id != '')
df$terL_coo[sele] <- sapply(df$terL_id[sele], function (x) {
                        i <- which(k$protein_id == x)
                        s <- ifelse(k$sign[i] == '+', 'f', 'r')
                        paste0(s, ';', k$start[i], ';', k$stop[i])
})



DF <- rbind(DF, df)




# ------------------------------ other databases -------------------------------------- #
pub <- c(
    MGV = 'Nayfach et al. 2021',
    GPD = 'Camarillo-Guerrero et al. 2021',
    GVD = 'Gregory et al. 2020',
    DEVoC = 'Van Espen et al. 2021',
    HuVirDB = 'Soto-Perez et al. 2019',
    BanfieldLab = ''
)



for (db in c('MGV', 'GPD', 'GVD', 'DEVoC', 'HuVirDB', 'BanfieldLab')) {

    t1 <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_TerL_genomes.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )



    t2 <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_SeleCode.txt'),
        sep = '\t',
        header = T,
        row.names = 1,
        stringsAsFactors = F
    )



    w <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG.txt'),
        header = T,
        sep = '\t',
        stringsAsFactors = F
    )



    k <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_proteomes.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )



    df <- data.frame(
        genome_id        = t1$genome_id,
        genome_len       = tab[t1$genome_id, 'Length'],
        genome_pgc       = tab[t1$genome_id, 'X.GC'],
        genome_code      = t2[t1$genome_id, 'sele_code'],
        prophage_cenote  = NA,
        genome_circ      = sapply(t1$genome_id, function (x) { term_rep_type(x, tab[x, 'Length'], w) }),
        cohort           = NA,
        source           = db,
        source_formal    = pub[db],
        terL_id          = t1$protein_id,
        terL_coo         = sapply(t1$protein_id, function (x) {
                                i <- which(k$protein_id == x)
                                s <- ifelse(k$strain[i] == 1, 'f', 'r')
                                paste0(s, ';', k$start[i], ';', k$end[i])
                           }),
        taxo_refseq      = '',
        taxo_Benler_2020 = '',
        taxo_BanfieldLab = '',
        taxo_our_crAss   = '',
        stringsAsFactors = F
    )



    DF <- rbind(DF, df)

}



v <- read.table('15_Lak_phage_genomes.ids', header = F, sep = ' ', stringsAsFactors = F)[, 1]
DF[DF$genome_id %in% v, 'source_formal'] <- 'Devoto et al. 2019'
DF[DF$genome_id %in% v, 'taxo_BanfieldLab'] <- 'Lak'



v <- read.table('349_huge_phage_genomes.ids', header = F, stringsAsFactors = F)[, 1]
DF[DF$genome_id %in% v, 'source_formal'] <- 'Al-Shayeb et al. 2020'



t <- read.table('/data/umcg-tifn/DATABASES/data_Al-Shayeb_2020/Table_S2_Clades.txt', sep = '\t', header = T, stringsAsFactors = F)

t <- t[!(t$Terminase.classification..group. %in% c('No terminase', 'Not classified', 'Not  classified', '')), ]

for (i in 1:nrow(t)) {

    x <- t$Genome.names.of.closely.related.genomes[i]
    y <- t$Terminase.classification..group.[i]

    if (x %in% DF$genome_id) { DF[DF$genome_id == x, 'taxo_BanfieldLab'] <- y }
}



t <- read.table('/data/umcg-tifn/DATABASES/data_Borges_2021_v1/borges_et_al_2021_phage_metadata.csv', sep = ',', header = T, stringsAsFactors = F)

DF[DF$genome_id %in% t$contig_id, 'source_formal'] <- 'Borges et al. 2021'

t <- t[t$Clade != 'No_clade_assigned', ]

for (i in 1:nrow(t)) {

    x <- t$contig_id[i]
    y <- t$Clade[i]

    if (x %in% DF$genome_id) { DF[DF$genome_id == x, 'taxo_BanfieldLab'] <- y }
}




# ------------------------------ species-level clusters ------------------------------ #
X <- read.table(
        'DB0_clusters_ANI95_AF85.tsv',
        sep = '\t',
        header = F,
        stringsAsFactors = F
)



DF$sp_repres <- ''

for (i in 1:nrow(X)) {

    v <- strsplit(X[i, 2], ',')[[1]]


    IDX <- which(DF$genome_id %in% v)
    idx <- IDX


    if (any(DF$genome_circ[idx] != '')) {    # priority to genomes with terminal repeats

        # select median length genome with terminal repeats
        idx <- idx[ DF$genome_circ[idx] != '' ]
        idx <- idx[ order(DF$genome_len[idx], decreasing = T) ]
        idx <- idx[ round(median(seq_along(idx))) ]
        DF$sp_repres[IDX] <- DF$genome_id[idx]

    } else {

        # priority to TerL-containing genome fragments
        if (any(DF$terL_id[idx] != '')) { idx <- idx[ DF$terL_id[idx] != '' ] }

        # select the longest genome (fragment) w/o terminal repeats
        idx <- idx[ order(DF$genome_len[idx], decreasing = T) ]
        idx <- idx[1]
        DF$sp_repres[IDX] <- DF$genome_id[idx]

    }

}




# ------------------------------ write table ------------------------------ #
write.table(
    DF,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB1_info.txt'
)
