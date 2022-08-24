sessionInfo()



ids <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', sep = '\t', header = F, stringsAsFactors = F)[, 1]



V <- c()

for (x in ids) {

    v <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot/', x, '_annot.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )$profile_name

    V <- c(V, v)

}

V <- sort(unique(V))



write.table(V, row.names = F, col.names = F, quote = F, file = 'db4_detected_doms.txt')
