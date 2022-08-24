sessionInfo()



# ------------------------------ viral refseq papilloma & polyoma recognized by relaxed cenote ------------------------------ #
DF <- read.table(
    'papilloma_polyoma_CONTIG_SUMMARY.tsv',
    sep = '\t',
    header = T,
    stringsAsFactors = F
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



# ------------------------------ % taxonomic group recognized by relaxed cenote ------------------------------ #
V <- c('Papillomaviridae', 'Polyomaviridae', 'Hepadnaviridae')

for (g in V) {

    n_all <- names(L)[unlist(lapply(L, function (v) g %in% v))]

    n_detected <- n_all[n_all %in% DF$ORIGINAL_NAME]

    p <- length(n_detected) / length(n_all) * 100
    p <- round(p, 2)

    cat('\n', p, '%', g, 'were recognized by relaxed Cenote-Taker2\n')

}
