sessionInfo()



sele <- read.table(
    'DB4_plus.ids',
    header = F,
    stringsAsFactors = F
)[, 1]

t0 <- read.table('/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt', sep = '\t', header = T, stringsAsFactors = F)



L <- list()



for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    cat('\nWorking with', k, 'cohort ...\n')


    samples <- t0$sample_id[ t0$cohort == k ]
    if (k == 'IBD') { samples <- read.table('/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/selected_ibd_samples.txt', header = F, stringsAsFactors = F)[, 1] }


    for (i in seq_along(samples)) {

        if (i %% 100 == 0) { cat('\tWorking with sample', i, '...\n') }
        s <- samples[i]


        f1 <- paste0('/data/umcg-tifn/NL_vir_analysis/map_reads/map_', k, '_reads/', s, '/', s, '.coverage.txt')
        if (!file.exists(f1)) { stop('File ', f1, ' is missing!') }
        t1 <- read.table(f1, sep = '\t', header = F, row.names = 1, stringsAsFactors = F)


        f2 <- paste0(s, '_DB4_cov_depth.txt')
        if (!file.exists(f2)) { stop('File ', f2, ' is missing!') }
        t2 <- read.table(f2, sep = '\t', header = F, stringsAsFactors = F)


        phages <- rownames(t1)[ t1[, 6] >= 0.75 ]
        phages <- phages[ phages %in% sele ]
        if (length(phages) == 0) { next }


        for (p in phages) {

            if (!(p %in% names(L))) { L[[ p ]] <- list() }
            if (!(k %in% names(L[[ p ]]))) { L[[ p ]][[ k ]] <- list() }


            len <- t1[p, 5]
            i <- which(t2[, 1] == p)
            coo <- t2[i, 2]
            cov <- t2[i, 3]


            COVERAGE <- rep(0, len)
            COVERAGE[coo] <- cov


            w_center <- seq(from = 501, to = len - 500, by = 200)
            w_from <- w_center - 500
            w_to <- w_center + 500
            w_to[ length(w_to) ] <- len


            L[[ p ]][[ k ]][[ s ]] <- sapply(seq_along(w_center), function (i) {

                mean( COVERAGE[ w_from[i]:w_to[i] ] )

            })

        }

    }

}



saveRDS(L, file = 'DB4_cov_depth.rds')
