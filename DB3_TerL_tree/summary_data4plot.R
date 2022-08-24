sessionInfo()




# -------------------------------------- read data -------------------------------------- #
ids3 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids3 <- rev(ids3)

ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids4 <- rev(ids4)

tab1 <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

tab2 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)




# -------------------------------------- genome length -------------------------------------- #
DF <- data.frame(
    sp_repres = ids3,
    repres_len = tab1[ids3, 'genome_len'],
    stringsAsFactors = F
)



for (db in c('DB3', 'DB4')) {

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }



    v <- setNames(tab1[ids, 'genome_len'], ids)
    v <- sort(v, decreasing = T)

    cat('\n\n\nLength of genomes representing vOTUs in', db, 'ranged from', min(v), 'to', max(v), 'nt, mean', round(mean(v)), 'nt\n\n')

    cat('The longest genome representing a vOTU in', db, 'was', names(v)[1], '\n\n')

    cat('The shortest genome representing a vOTU in', db, 'was', names(v)[length(v)], '\n\n')

    cat('There were', sum(v > 200000), 'genomes > 200 kb representing vOTUs in', db, '\n\n')



    v <- sapply(ids, function (x) {

        i <- which(tab1$sp_repres == x)

        l <- tab1[x, 'genome_len']

        b <- any((tab1$genome_circ[i] != '') & ((tab1$genome_len[i] < l*0.5) | (tab1$genome_len[i] > l*1.5)))

        return(b)

    })

    cat('Number of', db, 'vOTUs that contained genomes with terminal repeats that were > 50% smaller/larger than the representative genome:', sum(v), '\n')



    if (db == 'DB4') {

        Q <- data.frame(NULL)

        for (x in names(v)[v]) {

            i <- which(tab1$sp_repres == x)

            l <- tab1[x, 'genome_len']

            q <- data.frame(
                sp_id = x,
                n_len_same = sum((tab1$genome_circ[i] != '') & (tab1$genome_len[i] >= l*0.5) & (tab1$genome_len[i] <= l*1.5)),
                n_len_diff = sum((tab1$genome_circ[i] != '') & ((tab1$genome_len[i] < l*0.5) | (tab1$genome_len[i] > l*1.5))),
                stringsAsFactors = F
            )

            Q <- rbind(Q, q)

        }

        cat('\n'); print(Q, width = 200)

    }

}




# -------------------------------------- GC content -------------------------------------- #
DF$repres_pgc <- tab1[ids3, 'genome_pgc']



for (db in c('DB3', 'DB4')) {

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }

    v <- tab1[ids, 'genome_pgc']

    cat('\n\n\nGC content of genomes representing vOTUs in', db, 'ranged from', round(min(v)), 'to', round(max(v)), '%, mean', round(mean(v)), '%\n')

}




# -------------------------------------- terminal repeats -------------------------------------- #
DF$repres_TR <- tab1[ids3, 'genome_circ']



for (db in c('DB3', 'DB4')) {

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }



    v <- tab1[ids, 'genome_circ']

    cat('\n\n\nTermini of the genomes representing vOTUs in', db, ':\n')
    print(table(v, dnn = ''))
    cat('\n')



    v <- sapply(ids, function (x) {

        i <- which(tab1$sp_repres == x)

        r <- tab1[x, 'genome_circ']

        b <- any(!(tab1$genome_circ[i] %in% c('', r)))

        return(b)

    })

    cat('Number of', db, 'vOTUs that contained genomes with terminal repeats different from the representative genome:', sum(v), '\n')



    if (db == 'DB4') {

        Q <- data.frame(NULL)

        for (x in names(v)[v]) {

            i <- which(tab1$sp_repres == x)

            r <- tab1[x, 'genome_circ']

            q <- data.frame(
                sp_id = x,
                repres_TR = r,
                n_TR_same = sum(tab1$genome_circ[i] == r),
                n_TR_diff = sum(!(tab1$genome_circ[i] %in% c('', r))),
                stringsAsFactors = F
            )

            Q <- rbind(Q, q)

        }

        cat('\n'); print(Q, width = 200)

    }

}




# -------------------------------------- genetic code -------------------------------------- #
DF$repres_code <- tab1[ids3, 'genome_code']



for (db in c('DB3', 'DB4')) {

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }



    v <- tab1[ids, 'genome_code']

    cat('\n\n\nCode of the genomes representing vOTUs in', db, ':\n')
    print(table(v, dnn = ''))
    cat('\n')



    v <- sapply(ids, function (x) {

        i <- which(tab1$sp_repres == x)

        k <- tab1[x, 'genome_code']

        b <- any((tab1$genome_len[i] >= 20000) & (tab1$genome_code[i] != k))

        return(b)

    })

    cat('Number of', db, 'vOTUs that contained >= 20 kb genomes with code different from the representative genome:', sum(v), '\n\n')



    v <- sapply(ids, function (x) {

        i <- which(tab1$sp_repres == x)

        k <- tab1[x, 'genome_code']

        b <- any(tab1$genome_code[i] != k)

        return(b)

    })

    cat('Number of', db, 'vOTUs that contained genomes with code different from the representative genome:', sum(v), '\n')



    if (db == 'DB4') {

        Q <- data.frame(NULL)

        for (x in names(v)[v]) {

            i <- which(tab1$sp_repres == x)

            k <- tab1[x, 'genome_code']

            q <- data.frame(
                sp_id = x,
                repres_code = k,
                n_code_same = sum(tab1$genome_code[i] == k),
                n_code_diff = sum(tab1$genome_code[i] != k),
                n_code_same_from20kb = sum((tab1$genome_len[i] >= 20000) & (tab1$genome_code[i] == k)),
                n_code_diff_from20kb = sum((tab1$genome_len[i] >= 20000) & (tab1$genome_code[i] != k)),
                stringsAsFactors = F
            )

            Q <- rbind(Q, q)

        }

        cat('\n'); print(Q, width = 200)

    }

}




# -------------------------------------- prophages in vOTU -------------------------------------- #
DF$sp_prophage <- sapply(ids3, function (x) {

    i <- which(tab1$sp_repres == x)

    any((!is.na(tab1$prophage_cenote[i])) & (tab1$prophage_cenote[i] == 1)) * 1

})



for (db in c('DB3', 'DB4')) {

    v <- setNames(DF$sp_prophage, ids3)

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }

    cat('\n\n\n', sum(v[ids]), db, 'vOTUs contained prophage contigs\n')

}

cat('\n\n\n')




# -------------------------------------- % positive samples -------------------------------------- #
N <- c(LLD = 1135, LLD2 = 338, X300OB = 298, IBD = 455, DEVoC_all = 254, DEVoC_HA = 52)

for (k in names(N)) {

    k_ <- sub('^X', '', k)
    x1 <- paste0('n_samples_', k_)
    x2 <- paste0('p_samples_', k_)

    DF[, x2] <- tab2[ids3, x1] / N[k] * 100

    i <- which(DF[, x2] == max(DF[, x2]))[1]
    cat('The most prevalent phage in', k_, 'cohort was', DF$sp_repres[i], 'detected in', round(DF[i, x2]), '% samples\n\n')

}

cat('\n\n')



for (db in c('DB3', 'DB4')) {

    if (db == 'DB3') { ids <- ids3 } else { ids <- ids4 }

    n1 <- sum(tab2[ids, 'n_samples_DEVoC_all'] > 0)
    n2 <- sum(tab2[ids, 'n_samples_DEVoC_HA'] > 0)

    cat(n1, db, 'vOTUs were detected in DEVoC_all;', n2, db, 'vOTUs were detected in DEVoC_HA\n\n')

}

cat('\n\n')



write.table(DF, row.names = F, quote = F, sep = '\t', file = 'DB3_data4plot.txt')
