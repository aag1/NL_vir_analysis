.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(foreach)
library(meta)
sessionInfo()




# ------------------------------ read data ------------------------------ #
K <- c('LLD', 'LLD2', 'X300OB', 'IBD')


G <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]


L <- list()

for (k in K) {

    f <- paste0('/data/umcg-tifn/NL_vir_analysis/map_reads/', sub('^X', '', k), '_DB1_abundance.txt')
    L[[k]] <- read.table(f, row.names = 1, sep = '\t', header = T, stringsAsFactors = F)
    L[[k]] <- L[[k]][G, ]


    if (k == 'IBD') {

        S <- read.table('/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/selected_ibd_samples.txt', header = F, stringsAsFactors = F)[, 1]
        print(all(S %in% colnames(L[[k]])))
        L[[k]] <- L[[k]][, S]

    }

}


M <- read.table('/data/umcg-tifn/MetaPhlAn_4cohorts/LLD_LLD2_300OB_IBD_merged_abundance_table.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)
colnames(M) <- sub('_metaphlan$', '', colnames(M))




# ------------------------------ co-abundance analysis in individual cohorts ------------------------------ #
A <- list()

for (k in K) {

    # input data
    m1 <- L[[k]]
    m1 <- m1[apply(m1, 1, function (v) sum(v > 0)) > 10, ]

    m2 <- M[, colnames(m1)]
    m2 <- m2[apply(m2, 1, function (v) sum(v > 0)) > 10, ]

    # correlations
    DF <- foreach(i = 1:nrow(m1), .combine = rbind) %:%

        foreach(j = 1:nrow(m2), .combine = rbind) %do% {

            v1 <- as.numeric(m1[i, ])
            v2 <- as.numeric(m2[j, ])

            OBJ <- cor.test(x = v1, y = v2, method = 'spearman')

            data.frame(
                phage     = rownames(m1)[i],
                host_taxo = rownames(m2)[j],
                R         = OBJ$estimate,
                P         = OBJ$p.value,
                stringsAsFactors = F
            )

        }

    rownames(DF) <- paste0(DF$phage, ' ', DF$host_taxo)

    A[[k]] <- DF

}

warnings()




# ------------------------------ co-abundance meta-analysis ------------------------------ #
V1 <- unique(c(rownames(A[['LLD']]), rownames(A[['LLD2']]), rownames(A[['X300OB']]), rownames(A[['IBD']])))
V2 <- c('phage', 'host_taxo', 'LLD_R', 'LLD_P', 'LLD2_R', 'LLD2_P', 'X300OB_R', 'X300OB_P', 'IBD_R', 'IBD_P')

D <- matrix(NA, nrow = length(V1), ncol = length(V2), dimnames = list(V1, V2))
D <- as.data.frame(D)

for (k in K) {

    cR <- paste0(k, '_R')
    cP <- paste0(k, '_P')
    D[rownames(A[[k]]), c('phage', 'host_taxo', cR, cP)] <- A[[k]]

}


metas <- foreach(i = 1:nrow(D), .combine = rbind) %do% {

    Rs <- unlist(D[i, c('LLD_R', 'X300OB_R', 'IBD_R')])
    Ps <- unlist(D[i, c('LLD_P', 'X300OB_P', 'IBD_P')])
    sizes <- c(1135, 298, 455)

    sele <- which(!is.na(Rs))

    if (length(sele) == 0) {

        unlist(D[i, c('LLD2_R', 'LLD2_P')])

    } else if (length(sele) == 1) {

        c(Rs[sele], Ps[sele])

    } else { 

        OBJ <- metacor(Rs[sele], sizes[sele], sm = 'ZCOR', method.tau = 'SJ')
        c(OBJ$TE.random, OBJ$pval.random)

    }

}

warnings()


metas2 <- cbind(metas, p.adjust(metas[, 2], method = 'BH'))

colnames(metas2) <- c('meta_R', 'meta_P', 'meta_FDR')


D <- cbind(D, metas2)

D <- D[order(D$meta_FDR, decreasing = F), ]




# ------------------------------ write data ------------------------------ #
write.table(
    D,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB4_phage_host_coAbundance_metaAnalysis.txt'
)

idx <- sapply(G, function (g) which(D$phage == g)[1])

write.table(
    D[idx, ],
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB4_phage_host_coAbundance_metaAnalysis_top_corr.txt'
)
