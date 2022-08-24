.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
sessionInfo()




# ------------------------------ read data ------------------------------ #
t1 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)


t2 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)


tr <- read.tree('/data/umcg-tifn/NL_vir_analysis/DB3_TerL_tree/DB3_terL_tree.rooted.newick')




# ------------------------------ DB3 ids ------------------------------ #
ids3 <- rev(tr$tip.label)

write.table(ids3, row.names = F, col.names = F, quote = F, file = 'DB3.tree.ids')




# ------------------------------ DB4 ids ------------------------------ #
t <- t2[ids3, ]

t <- t[(t$n_samples_LLD / 1135 > 0.05) | (t$n_samples_LLD2 / 338 > 0.05) | (t$n_samples_300OB / 298 > 0.05) | (t$n_samples_IBD / 455 > 0.05), ]

ids4 <- rownames(t)

write.table(ids4, row.names = F, col.names = F, quote = F, file = 'DB4.tree.ids')




# ------------------------------ DB3 table 1 ------------------------------ #
IDX <- c()

for (x in ids3) {

    idx <- which(t1$sp_repres == x)
    idx <- idx[ order(t1$genome_len[idx], decreasing = T) ]
    IDX <- c(IDX, idx)

}

t <- t1[IDX, ]


d <- t[, c('sp_repres', 'genome_id', 'genome_len', 'genome_pgc', 'genome_code', 'genome_circ', 'prophage_cenote', 'source_formal', 'taxo_our_crAss', 'taxo_Benler_2020', 'taxo_refseq')]


d$genome_code <- sub('^c', '', d$genome_code)


d$prophage_cenote[ is.na(d$prophage_cenote) ] <- ''
d$prophage_cenote[ d$prophage_cenote == 0 ] <- ''
d$prophage_cenote[ d$prophage_cenote == 1 ] <- 'yes'


d$taxo_Benler_2020 <- sapply(d$taxo_Benler_2020, function (x) {

    v <- strsplit(x, ',')[[1]]
    f2 <- grep('viridae$', v, value = T)
    if (length(f2) == 0) { f2 <- '' } else { f2 <- f2[1] }
    return(f2)

})


write.table(d, sep = '\t', row.names = F, quote = F, file = 'DB3_table1.txt')




# ------------------------------ NL_vir from DB3 vOTUs ------------------------------ #
ids <- t$genome_id[ t$source == 'NL_vir' ]

write.table(ids, row.names = F, col.names = F, quote = F, file = 'NL_vir_from_DB3_vOTUs.ids')




# ------------------------------ DB3 table 2 ------------------------------ #
cat('\n\n')
t1[(t1$genome_id %in% ids3) & (t1$source == 'RefSeq'), c('genome_id', 'terL_id', 'terL_coo')]
cat('\n\n')


q <- data.frame(
    sp_repres = ids3,
    sp_repres_terL = sapply(ids3, function (x) {
                                        i <- which(t1$genome_id == x)
                                        ifelse(
                                            t1$source[i] == 'RefSeq',
                                            sub('^[^:]+:([0-9]+)\\.\\.([0-9]+)$', 'f;\\1;\\2', t1$terL_coo[i]),
                                            t1$terL_coo[i]
                                        )}),
    stringsAsFactors = F
)


q <- cbind(q, t2[ids3, 1:6])


write.table(q, sep = '\t', row.names = F, quote = F, file = 'DB3_table2.txt')
