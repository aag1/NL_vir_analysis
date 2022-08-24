sessionInfo()



G <- unique(read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', header = T, stringsAsFactors = F)$sp_repres)



for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    f1 <- paste0('/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/', k, '_raw_clean_reads_number.txt')
    t1 <- read.table(f1, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
    S <- rownames(t1)
    if (k == 'IBD') { S <- read.table('/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/selected_ibd_samples.txt', header = F, stringsAsFactors = F)[, 1] }


    DF <- data.frame(NULL)
    M <- matrix(0, nrow = length(G), ncol = length(S), dimnames = list(G, S))


    for (x in S) {

        f2 <- paste0('/data/umcg-tifn/NL_vir_analysis/map_reads/map_', k, '_reads/', x, '/', x, '.coverage.txt')
        t2 <- read.table(f2, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
        t2 <- t2[G, ]
        idx <- which(t2[, 6] >= 0.75)


        df <- data.frame(
            sample_id = x,
            num_sp_cov75 = length(idx),
            pct_reads_mapped = round(sum(t2[, 3]) / t1[x, 'clean_reads'] * 100, 3),
            pct_reads_mapped_cov75 = round(sum(t2[idx, 3]) / t1[x, 'clean_reads'] * 100, 3),
            stringsAsFactors = F
        )
        DF <- rbind(DF, df)


        if (length(idx) == 0) { next }
        M[idx, x] <- t2[idx, 3] / t2[idx, 5] / t1[x, 'clean_reads'] * 10^6

    }


    write.table(DF, sep = '\t', quote = F, row.names = F, file = paste0(k, '_DB1_mapping_summary.txt'))
    write.table(M, sep = '\t', quote = F, file = paste0(k, '_DB1_abundance.txt'))
    
}
