.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()




# ------------------------------ input data ------------------------------ #
ids <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', sep = '\t', header = F, stringsAsFactors = F)[, 1]



### proteome coo & seq
G <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

S <- read.fasta(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.fasta',
    seqtype = 'AA',
    as.string = T,
    set.attributes = F
)


for (db in c('MGV', 'GPD', 'GVD', 'DEVoC', 'HuVirDB', 'BanfieldLab')) {

    t <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_proteomes.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )

    G <- rbind(G, t)

    x <- read.fasta(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_proteomes.fasta'),
        seqtype = 'AA',
        as.string = T,
        set.attributes = F
    )

    S <- c(S, x)

}


t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

G <- rbind(G, t)

x <- read.fasta(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.fasta',
    seqtype = 'AA',
    as.string = T,
    set.attributes = F
)

S <- c(S, x)



### proteome annot
DF <- data.frame(NULL)

for (x in ids) {

    df <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot/', x, '_annot.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}




# ------------------------------ extract selected proteins ------------------------------ #
for (x in c('RVT_1', 'Phage_integrase', 'rve')) {

    p_id <- DF$protein_id[ DF$profile_name == x ]


    tab <- G[sapply(p_id, function (x) which(G$protein_id == x)), ]

    tab$orf_coo <- paste0(ifelse(tab$strain == 1, 'f', 'r'), ';', tab$start, ';', tab$end)

    tab$genome_id_short <- sapply(tab$genome_id, function (x) sub('^(.+_NODE_[0-9]+)_length_[0-9]+_cov_[0-9\\.]+$', '\\1', x))

    tab$protein_id_new <- paste0(tab$genome_id_short, '_', tab$orf_coo)

    tab <- tab[, c('genome_id', 'protein_id', 'protein_id_new', 'orf_coo')]

    write.table(tab, sep = '\t', quote = F, row.names = F, file = paste0('DB4_', x, '_info.txt'))


    seq <- S[ tab$protein_id ]

    seq <- lapply(seq, function (x) sub('\\*$', '', x))

    write.fasta(seq, names = tab$protein_id_new, file.out = paste0('DB4_', x, '_seq.fasta'), nbchar = 80, as.string = T)

}
