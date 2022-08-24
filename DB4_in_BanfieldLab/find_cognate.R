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
                db_id    = x,
                db_len   = tab$slen[idx][1],
                db_cov   = round(sCOV * 100, 2),
                db_rev   = (tab$sstrand[idx][1] == 'minus') * 1,
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

for (x in V) {

    for (y in c('Devoto_2019', 'Al-Shayeb_2020', 'Borges_2021')) {

        f <- paste0(x, '_vs_', y, '.txt')
        if (file.info(f)$size == 0) { next }

        t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)
        colnames(t) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')

        t1 <- find_cognate(x, t, 0.5, 0.1, 0.5)
        if (nrow(t1) != 0) { t1$db_name <- y }

        if (nrow(t1) != 0) { if (nrow(T1) == 0) { T1 <- t1 } else { T1 <- rbind(T1, t1) } }

    }

}




# ------------------------------ cognate info ------------------------------ #
d <- read.table('/data/umcg-tifn/DATABASES/data_Borges_2021_v1/borges_et_al_2021_phage_metadata.csv', sep = ',', header = T, stringsAsFactors = F)

T1$Borges_2021_info <- NA

sele <- which(T1$db_id %in% d$contig_id)

T1$Borges_2021_info[sele] <- sapply(T1$db_id[sele], function (x) {

    paste(d[d$contig_id == x, c('Clade', 'Group', 'ecosystem', 'original_scaffold_prophage_status', 'gcode')], collapse = ',')

})




# ------------------------------ write table ------------------------------ #
write.table(T1, sep = '\t', quote = F, row.names = F, file = 'DB4_cognate_BanfieldLab.txt')
