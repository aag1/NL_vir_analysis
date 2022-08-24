sessionInfo()




# ------------------------------ read data ------------------------------ #
tab <- read.table('DB3_crispr_links.txt', sep = '\t', header = T, stringsAsFactors = F)

d <- read.table('DB3_crispr_taxo.txt', sep = '\t', header = F, stringsAsFactors = F)
colnames(d) <- c('host', 'host_taxo', 'host_genus')

V <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]




# ------------------------------ assign taxonomy ------------------------------ #
tab[, c('host_taxo', 'host_genus')] <- ''

for (i in 1:nrow(d)) {

    idx <- which(tab$host == d$host[i])
    tab[idx, 'host_taxo'] <- d$host_taxo[i]
    tab[idx, 'host_genus'] <- d$host_genus[i]

}




# ------------------------------ write data ------------------------------ #
tab <- tab[, c('phage', 'host_taxo', 'host_genus')]

tab <- unique(tab)

IDX <- c()

for (x in V) {

    idx <- which(tab$phage == x)

    IDX <- c(IDX, idx)

}

write.table(tab[IDX, ], sep = '\t', quote = F, row.names = F, file = 'DB3_crispr_hosts.txt')
