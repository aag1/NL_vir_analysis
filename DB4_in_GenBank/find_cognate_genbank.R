.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(IRanges)
sessionInfo()




# ------------------------------ function ------------------------------ #
find_cognate <- function(query, tab, ident_thr = 0.95, qcov_thr = 0.1, scov_thr = 0.85) {

    tab <- tab[tab$pident >= ident_thr, ]

    DF <- data.frame(NULL)

    for (x in unique(tab$sseqid)) {

        idx <- which(tab$sseqid == x)

        qFROM <- c()
        qTO <- c()

        sFROM <- c()
        sTO <- c()

        for (i in idx) {

            qFROM <- c(qFROM, tab$qstart[i])
            qTO   <- c(qTO,   tab$qend[i])

            sFROM <- c(sFROM, ifelse(tab$sstrand[i] == 'plus', tab$sstart[i], tab$send[i]))
            sTO   <- c(sTO,   ifelse(tab$sstrand[i] == 'plus', tab$send[i], tab$sstart[i]))

        }

        qIR  <- IRanges(start = qFROM, end = qTO)
        qCL  <- reduce(qIR)
        qCL  <- as.data.frame(qCL)
        qCOV <- sum(qCL$width) / tab$qlen[idx][1]

        if (qCOV < qcov_thr) { next }

        sIR  <- IRanges(start = sFROM, end = sTO)
        sCL  <- reduce(sIR)
        sCL  <- as.data.frame(sCL)
        sCOV <- sum(sCL$width) / tab$slen[idx][1]

        if (sCOV >= scov_thr) {

            df <- data.frame(
                nl_id    = query,
                nl_len   = tab$qlen[idx][1],
                nl_cov   = round(qCOV * 100, 2),
                gb_id    = strsplit(x, '\\|')[[1]][4],
                gb_len   = tab$slen[idx][1],
                gb_cov   = round(sCOV * 100, 2),
                gb_name  = tab$sscinames[idx][1],
                gb_title = tab$stitle[idx][1],
                gb_rev   = (tab$sstrand[idx][1] == 'minus') * 1,
                stringsAsFactors = F
            )

            df <- df[order(df$nl_cov, decreasing = T), ]

            if (nrow(DF) == 0) { DF <- df } else { DF <- rbind(DF, df) }

        }

    }

    return(DF)

}




# ------------------------------ identify cognate ------------------------------ #
V <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]

T1 <- data.frame(NULL)
T2 <- data.frame(NULL)

for (x in V) {

    f <- paste0(x, '_vs_GenBank.txt')

    t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)
    colnames(t) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident', 'staxids', 'sscinames', 'stitle')

    t1 <- find_cognate(x, t, 0.95, 0.1, 0.85)
    t2 <- find_cognate(x, t, 0.5, 0.1, 0.5)

    if (nrow(t1) != 0) { if (nrow(T1) == 0) { T1 <- t1 } else { T1 <- rbind(T1, t1) } }
    if (nrow(t2) != 0) { if (nrow(T2) == 0) { T2 <- t2 } else { T2 <- rbind(T2, t2) } }

}




# ------------------------------ write tables ------------------------------ #
write.table(T1, sep = '\t', quote = F, row.names = F, file = 'DB4_cognate_GenBank_species_level.txt')

write.table(T2, sep = '\t', quote = F, row.names = F, file = 'DB4_cognate_GenBank_genera_level.txt')
