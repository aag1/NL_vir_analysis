sessionInfo()




# ------------------------------ number of clean reads per sample ------------------------------ #
f1 <- 'DEVoC_clean_reads_multiqc_report_data/mqc_fastqc_sequence_counts_plot_1.txt'
t1 <- read.table(f1, sep = '\t', header = T, stringsAsFactors = F)


tab <- data.frame(
    sample_id = unique(sub('^([^_]+)_kneaddata_paired_[12]$', '\\1', t1$Sample)),
    clean_reads = NA,
    stringsAsFactors = F
)
rownames(tab) <- tab$sample_id


for (x in tab$sample_id) {

    i1 <- which(t1$Sample == paste0(x, '_kneaddata_paired_1'))
    i2 <- which(t1$Sample == paste0(x, '_kneaddata_paired_2'))

    n1 <- sum(t1$Unique.Reads[i1] + t1$Duplicate.Reads[i1])
    n2 <- sum(t1$Unique.Reads[i2] + t1$Duplicate.Reads[i2])
    if (n1 != n2) { stop('n1 != n2 for ', x) }

    tab[x, 'clean_reads'] <- n1 + n2

}


S <- tab$sample_id




# ------------------------------ genome IDs ------------------------------ #
G <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.ids', header = F, stringsAsFactors = F)[, 1]




# ------------------------------ coverage ------------------------------ #
DF <- data.frame(NULL)


for (x in S) {

    f2 <- paste0(x, '.coverage.txt')
    t2 <- read.table(f2, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
    t2 <- t2[G, ]
    idx <- which(t2[, 6] >= 0.75)


    df <- data.frame(
        sample_id = x,
        num_sp_cov75 = length(idx),
        pct_reads_mapped = round(sum(t2[, 3]) / tab[x, 'clean_reads'] * 100, 3),
        pct_reads_mapped_cov75 = round(sum(t2[idx, 3]) / tab[x, 'clean_reads'] * 100, 3),
        stringsAsFactors = F
    )
    DF <- rbind(DF, df)

}


write.table(DF, sep = '\t', quote = F, row.names = F, file = 'DEVoC_DB3_mapping_summary.txt')
