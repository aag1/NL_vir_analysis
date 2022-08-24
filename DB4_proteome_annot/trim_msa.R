.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
sessionInfo()



for (x in c('RVT_1', 'Phage_integrase', 'rve')) {

    t <- read.table(paste0('DB4_', x, '_info.txt'), sep = '\t', header = T, stringsAsFactors = F)

    a <- read.fasta(paste0(x, '_msa_with_pfam.fasta'))$ali

    a <- a[t$protein_id_new, ]

    empty_cols <- apply(a, 2, function (v) all(v == '-'))

    a <- a[, !empty_cols]

    write.fasta(ids = rownames(a), seqs = a, file = paste0('DB4_', x, '_msa.fasta'))

}
