.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(dplyr)
library(alluvial)
sessionInfo()




# ------------------------------ input data ------------------------------ #
DF <- read.table(
    'DB2_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)

DF$simple_termini <- ifelse(DF$simple_termini == 'TR', 'Yes', 'No')

DF$simple_taxo[ !(DF$simple_taxo %in% c('Unassigned', 'Caudoviricetes')) ] <- 'Other'




# ------------------------------ sankey plot ------------------------------ #
pdf('DB2_sankey.pdf', height = 3.35, width = 3.35)


P <- setNames(c('blue', 'red', 'grey80'), c('Caudoviricetes', 'Other', 'Unassigned'))


A <- DF[, c('simple_taxo', 'simple_termini', 'simple_source')] %>% group_by_all() %>% summarise(freq = n())
A <- as.data.frame(A)


alluvial(
    A[, c('simple_taxo', 'simple_termini', 'simple_source')],
    freq = A$freq,
    axis_labels = c('Taxonomy', 'Terminal repeats', 'Source'),
    col = P[ A$simple_taxo ],
    alpha = ifelse(A$simple_taxo == 'Other', 1, 0.5),
    border = sapply(A$simple_taxo, function (x) adjustcolor(P[x], alpha.f = ifelse(x == 'Other', 1, 0.5))),
    cw = 0.25,
    cex = 0.5, cex.axis = 0.5
)

dev.off()
