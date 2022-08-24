.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(IRanges)
sessionInfo()



tab <- read.table('crAss_sensu_stricto_vs_MGV-GENOME-0359371.txt', sep = '\t', header = T, stringsAsFactors = F)
colnames(tab) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')



ir <- IRanges(start = tab$qstart, end = tab$qend)
cl <- reduce(ir)
cl <- as.data.frame(cl)



n <- sum(cl$width) / tab$qlen[1] * 100
cat('\n\n\nQuery length covered by hits:', round(n), '%\n\n\n')
