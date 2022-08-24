sessionInfo()




# ------------------------------ read data ------------------------------ #
t1 <- read.table(
    'DB3_host_fragm_hosts.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)


t2 <- read.table(
    'DB3_crispr_hosts.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)


t3 <- read.table(
    'DB4_phage_host_coAbundance_metaAnalysis_top_corr.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ build table ------------------------------ #
DF <- data.frame(NULL)


for (x in t3$phage) {

    v1 <- t1$host_taxo[ t1$sp_id == x ]
    v2 <- t2$host_taxo[ t2$phage == x ]
    v3 <- t3$host_taxo[ t3$phage == x ]

    if ((length(v1) == 0) & (length(v2) == 0)) { next }

    v3 <- sub('^[a-z]__', '', v3)
    v3 <- gsub('\\|[a-z]__', ', ', v3)

    df <- data.frame(
        phage = x,
        host_source = c(rep('prophage', length(v1)), rep('CRISPR', length(v2)), 'co-abundance'),
        host_taxo = c(v1, v2, v3),
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}


write.table(
    DF,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB4_compare_host_predictions.txt'
)
