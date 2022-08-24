.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(IRanges)
sessionInfo()




# ------------------------------ collect info, select hits >= 1000 nt ------------------------------ #
DF <- data.frame(NULL)

V <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident', 'staxids', 'sscinames')

files <- list.files(pattern = '_host_fragm_vs_ncbi.txt$')

for (f in files) {

    if (file.info(f)$size == 0) { next }
    df <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)
    colnames(df) <- V

    df$sp_id <- sub('_host_fragm_vs_ncbi.txt$', '', f)
    df <- df[, c('sp_id', V)]

    df$sseqid <- unlist(lapply(strsplit(df$sseqid, '\\|'), function (v) v[4]))

    df <- df[df$length >= 1000, ]

    DF <- rbind(DF, df)

}




# ------------------------------ per query, select target producing maximal query coverage ------------------------------ #
IDX <- c()

for (x in unique(DF$qseqid)) {

    idx1 <- which(DF$qseqid == x)

    Q <- data.frame(NULL)

    for (y in unique(DF$sseqid[idx1])) {

        idx2 <- idx1[ DF$sseqid[idx1] == y ]

        ir <- IRanges(start = DF$qstart[idx2], end = DF$qend[idx2])
        cl <- reduce(ir)
        cl <- as.data.frame(cl)

        q <- data.frame(sseqid = y, qcov = sum(cl$width), stringsAsFactors = F)
        Q <- rbind(Q, q)

    }

    sele <- Q$sseqid[ Q$qcov == max(Q$qcov) ][1]
    idx3 <- idx1[ DF$sseqid[idx1] == sele ]
    IDX <- c(IDX, idx3)

}

DF <- DF[IDX, ]

write.table(DF, sep = '\t', quote = F, row.names = F, file = 'DB3_host_fragm_sele_hits.txt')
