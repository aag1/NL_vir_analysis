sessionInfo()



# ------------------------------ viral refseq genomes recognized by cenote ------------------------------ #
DF <- data.frame(NULL)

for (i in 0:99) {

    f <- paste0('viral_refseq_209_', i, '_vir_CONTIG_SUMMARY.tsv')

    df <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
    DF <- rbind(DF, df)

}

write.table(
    DF,
    sep = '\t',
    row.names = F,
    quote = F,
    'vir_refseq_recognized_by_cenote.txt'
)



# ------------------------------ viral refseq taxonomy ------------------------------ #
tab <- read.table(
    'viral_refseq_taxo.txt',
    sep = '\t',
    header = T,
    quote = '',
    fill = T,
    stringsAsFactors = F
)

L <- sapply(tab$genome_id, function (x) {

    v <- strsplit(tab[tab$genome_id == x, 'taxonomy'], ';')[[1]]
    if (length(v) == 0) { v <- '' }
    return(v)

})



# ------------------------------ % taxonomic group recognized by cenote ------------------------------ #
V <- c(
    'Caudoviricetes',
    'Herpesviridae',
    'Papillomaviridae',
    'Polyomaviridae',
    'Adenoviridae',
    'Tectiliviricetes_without_Adenoviridae',
    'Nucleocytoviricota',
    'Hepadnaviridae'
)

for (g in V) {

    if (g == 'Tectiliviricetes_without_Adenoviridae') {

        n_all <- names(L)[unlist(lapply(L, function (v) ('Tectiliviricetes' %in% v) & (!('Adenoviridae' %in% v))))]

    } else {

        n_all <- names(L)[unlist(lapply(L, function (v) g %in% v))]

    }


    n_detected <- n_all[n_all %in% DF$ORIGINAL_NAME]


    p <- length(n_detected) / length(n_all) * 100
    p <- round(p, 2)


    cat('\n', p, '%', g, 'were recognized by Cenote-Taker2\n')

}
