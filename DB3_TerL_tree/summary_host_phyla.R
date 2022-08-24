sessionInfo()




# ------------------------------ read data ------------------------------ #
ids3 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids3 <- rev(ids3)

ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids4 <- rev(ids4)

tab1 <- read.table('/data/umcg-tifn/NL_vir_analysis/host_prediction/DB3_host_fragm_hosts.txt', sep = '\t', header = T, stringsAsFactors = F)

tab2 <- read.table('/data/umcg-tifn/NL_vir_analysis/host_prediction/DB3_crispr_hosts.txt', sep = '\t', header = T, stringsAsFactors = F)




# ------------------------------ create data frame ------------------------------ #
DF <- data.frame(
    host_phylum_prophage = rep('', length(ids3)),
    host_phylum_crispr = '',
    host_phylum_consensus = '',
    stringsAsFactors = F
)

rownames(DF) <- ids3




# ------------------------------ host phyla based on prophages ------------------------------ #
tab1$host_phylum <- unlist(lapply(strsplit(tab1$host_taxo, ','), function (v) {

    v <- sub('^ ', '', v)

    ifelse(length(v) < 2, '', v[2])

}))

tab1$host_phylum[tab1$host_phylum %in% c('bacterium', 'environmental samples', 'environmental samples Bacteria')] <- ''



n <- 0

for (x in unique(tab1$sp_id)) {

    p <- unique(tab1$host_phylum[ tab1$sp_id == x ])

    if (length(p) > 1) { p <- p[p != ''] }

    if (length(p) > 1) { n <- n + 1; next }

    DF[x, 'host_phylum_prophage'] <- p

}



cat('\n\n\nProphage-based host predictions were available for', sum(DF$host_phylum_prophage != ''), 'vOTUs.\n\n')

cat('There were', n, 'vOTUs linked to multiple host phyla based on prophages; prophage-based host predictions for these vOTUs were disregarded.\n\n\n\n')




# ------------------------------ host phyla based on CRISPR ------------------------------ #
tab2$host_phylum <- unlist(lapply(strsplit(tab2$host_taxo, ','), function (v) {

    v <- sub('^ ', '', v)

    ifelse(length(v) < 2, '', v[2])

}))

tab2$host_phylum[tab2$host_phylum == 'bacterium LF-3'] <- ''



n <- 0

for (x in unique(tab2$phage)) {

    p <- unique(tab2$host_phylum[ tab2$phage == x ])

    if (length(p) > 1) { p <- p[p != ''] }

    if (length(p) > 1) { n <- n + 1; next }

    DF[x, 'host_phylum_crispr'] <- p

}



cat('CRISPR-based host predictions were available for', sum(DF$host_phylum_crispr != ''), 'vOTUs.\n\n')

cat('There were', n, 'vOTUs linked to multiple host phyla based on CRISPR; CRISPR-based host predictions for these vOTUs were disregarded.\n\n\n\n')




# ------------------------------ consensus host phyla ------------------------------ #
idx <- which((DF$host_phylum_prophage != '') & (DF$host_phylum_crispr == ''))

DF$host_phylum_consensus[idx] <- DF$host_phylum_prophage[idx]



idx <- which((DF$host_phylum_prophage == '') & (DF$host_phylum_crispr != ''))

DF$host_phylum_consensus[idx] <- DF$host_phylum_crispr[idx]



idx <- which((DF$host_phylum_prophage != '') & (DF$host_phylum_crispr != ''))

DF$host_phylum_consensus[idx] <- DF$host_phylum_prophage[idx]



write.table(DF, quote = F, sep = '\t', file = 'DB3_predicted_host_phyla.txt')



cat('Both prophage-based and CRISPR-based host predictions were available for', length(idx), 'vOTUs.\n\n')

bad <- idx[ DF$host_phylum_prophage[idx] != DF$host_phylum_crispr[idx] ]
cat('Conflicts between prophage-based and CRISPR-based host predictions:\n\n')
DF[bad, 1:2]
cat('\n\n\n')



cat('In total, host predictions were available for', sum(DF$host_phylum_consensus != ''), 'DB3 vOTUs:\n')
sort(table(DF$host_phylum_consensus), decreasing = T)
cat('\n\n\n')



DF <- DF[ids4, ]
cat('In total, host predictions were available for', sum(DF$host_phylum_consensus != ''), 'DB4 vOTUs:\n')
sort(table(DF$host_phylum_consensus), decreasing = T)
cat('\n\n\n')
