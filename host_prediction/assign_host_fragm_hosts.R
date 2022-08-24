sessionInfo()




# ------------------------------ read data ------------------------------ #
tab <- read.table('DB3_host_fragm_sele_hits.txt', sep = '\t', header = T, stringsAsFactors = F)

d <- read.table('DB3_host_fragm_taxo.txt', sep = '\t', header = F, stringsAsFactors = F)
colnames(d) <- c('sseqid', 'host_taxo', 'host_genus')

V <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]




# ------------------------------ assign taxonomy ------------------------------ #
tab[, c('host_taxo', 'host_genus')] <- ''

for (i in 1:nrow(d)) {

    idx <- which(tab$sseqid == d$sseqid[i])
    tab[idx, 'host_taxo'] <- d$host_taxo[i]
    tab[idx, 'host_genus'] <- d$host_genus[i]

}




# ------------------------------ write data ------------------------------ #
tab <- tab[, c('sp_id', 'host_taxo', 'host_genus')]

tab <- unique(tab)

IDX <- c()

for (x in V) {

    idx <- which(tab$sp_id == x)
    
    IDX <- c(IDX, idx)

}

write.table(tab[IDX, ], sep = '\t', quote = F, row.names = F, file = 'DB3_host_fragm_hosts.txt')
