sessionInfo()


d <- read.table('DB4_cognate_GenBank_pub_info.txt', sep = '\t', header = T, stringsAsFactors = F)


for (x in c('species', 'genera')) {

    f <- paste0('DB4_cognate_GenBank_', x, '_level.txt')
    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)

    q <- merge(t, d, by = 'gb_id', all.x = T, sort = F)
    q <- q[, c(2:4, 1, 5:8, 10:11, 9)]

    f <- paste0('DB4_cognate_GenBank_', x, '_level_with_pub_info.txt')
    write.table(q, sep = '\t', quote = F, row.names = F, file = f)

}
