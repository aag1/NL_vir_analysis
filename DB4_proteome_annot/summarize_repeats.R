sessionInfo()



files <- list.files(path = '/data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot/', pattern = '_repeats.6.txt', full.names = T)

DF <- data.frame(NULL)

for (f in files) {

    t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)
    colnames (t) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')


    bad <- c()
    for (i in 1:nrow(t)) {

        # exclude the hit between the whole genome and itself
        if (t$length[i] == t$qlen[1]) { bad <- c(bad, i) }

        # exclude terminal repeats
        if ((min(t$qstart[i], t$qend[i]) == 1) & (max(t$sstart[i], t$send[i]) == t$qlen[1])) { bad <- c(bad, i) }
        if ((max(t$qstart[i], t$qend[i]) == t$qlen[1]) & (min(t$sstart[i], t$send[i]) == 1)) { bad <- c(bad, i) }

        # exclude duplicates
        if (i > 1) {

            for (j in 1:(i-1)) {

                if (t$sstrand[i] != t$sstrand[j]) { next }

                if ((t$sstrand[i] == 'plus') &
                    (t$qstart[i] == t$sstart[j]) & (t$qend[i] == t$send[j]) &
                    (t$sstart[i] == t$qstart[j]) & (t$send[i] == t$qend[j]) ) { bad <- c(bad, i) }

                if ((t$sstrand[i] == 'minus') &
                    (t$qstart[i] == t$send[j]) & (t$qend[i] == t$sstart[j]) &
                    (t$send[i] == t$qstart[j]) & (t$sstart[i] == t$qend[j]) ) { bad <- c(bad, i) }

            }

        }

        # exclude short alignments & low identity alignments
        if ((t$length[i] < 100) | (t$pident[i] < 80)) { bad <- c(bad, i) }

    }


    if (t$qseqid[1] == 'ERS698742_NODE_4_length_57489_cov_840.739161') { bad <- c(bad, 1:2) }
    if (length(bad) > 0) { t <- t[-bad, ] }


    if (nrow(t) > 0) {
        cat(nrow(t), 'pair(s) of repeats detected in', t$qseqid[1], '\n\n\n')
        DF <- rbind(DF, t)
    }

}



DF$coo1 <- sapply(1:nrow(DF), function (i) (DF$qstart[i] + DF$qend[i]) / 2)
DF$coo2 <- sapply(1:nrow(DF), function (i) (DF$sstart[i] + DF$send[i]) / 2)

DF$from <- sapply(1:nrow(DF), function (i) min(DF$coo1[i], DF$coo2[i]))
DF$to   <- sapply(1:nrow(DF), function (i) max(DF$coo1[i], DF$coo2[i]))

DF <- DF[order(DF$qseqid, DF$from, DF$to), ]



write.table(DF, sep = '\t', quote = F, row.names = F, file = 'DB4_repeats.txt')
