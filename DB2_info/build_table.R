sessionInfo()




# ------------------------------ read data ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)


DF <- data.frame(NULL)

for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    f <- paste0('/data/umcg-tifn/NL_vir_analysis/map_reads/', k, '_DB1_abundance.txt')
    d <- read.table(f, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

    df <- data.frame(apply(d, 1, function (v) sum(v > 0)), stringsAsFactors = F)
    colnames(df) <- paste0('n_samples_', k)

    if (ncol(DF) == 0) { DF <- df } else { DF <- cbind(DF, df) }

}


f <- '/data/umcg-tifn/NL_vir_analysis/DEVoC_reads/DEVoC_DB1_abundance.txt'
d <- read.table(f, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

f <- '/home/umcg-agulyaeva/NL_vir_analysis/DEVoC_reads/DEVoC_adult_cohort_PRJNA722819.tsv'
q <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
h <- q$run_accession[ grepl('^HA', q$sample_alias) ]

df <- data.frame(
    n_samples_DEVoC_all = apply(d, 1, function (v) sum(v > 0)),
    n_samples_DEVoC_HA = apply(d[, h], 1, function (v) sum(v > 0)),
    stringsAsFactors = F
)

DF <- cbind(DF, df)

DF <- round(DF, 2)




# ------------------------------ select vOTUs detected in at least one Dutch sample ------------------------------ #
sele <- rownames(DF)[ apply(DF, 1, function (v) any(v[1:4] > 0)) ]

n1 <- sum(t$source[ t$sp_repres %in% sele ] == 'NL_vir')
n2 <- sum(t$source[ !(t$sp_repres %in% sele) ] == 'NL_vir')

cat('\nTotal NL_vir contigs:', n1 + n2, '\n')
cat('\nNL_vir contigs belonging to detected vOTUs:', n1, '(', round(n1/(n1+n2)*100, 2), '% )\n')
cat('\nNL_vir contigs belonging to undetected vOTUs:', n2, '(', round(n2/(n1+n2)*100, 2), '% )\n')

DF <- DF[sele, ]




# ------------------------------ simple taxo ------------------------------ #
DF$simple_taxo <- sapply(rownames(DF), function (x) {

    idx <- which(t$sp_repres == x)


    taxo <- c()


    ### refseq taxo
    W <- t$taxo_refseq[idx]
    W <- unique(W)
    W <- W[ W != '' ]

    if (length(W) > 0) {

        for (y in W) {

            v <- strsplit(y, ';')[[1]]

            if (any(v == 'Caudoviricetes')) {

                taxo <- c(taxo, 'Caudoviricetes')

            } else {

                if (any(grepl('viridae$', v))) {

                    taxo <- c(taxo, grep('viridae$', v, value = T)[1])
    
                } else {

                    taxo <- c(taxo, y)

                }
            }
        }
    }


    ### databases caudo
    if (any(t$source[idx] %in% c('Benler_2020', 'MGV', 'GPD', 'GVD', 'DEVoC', 'HuVirDB', 'BanfieldLab'))) {

        taxo <- c(taxo, 'Caudoviricetes')

    }


    ### NL_vir caudo
    if (any((t$source[idx] == 'NL_vir') & (t$terL_id[idx] != ''))) {

        taxo <- c(taxo, 'Caudoviricetes')

    }


    ### NL_vir adeno
    if ('NL_vir016001' %in% t$genome_id[idx]) {

        taxo <- c(taxo, 'Adenoviridae')

    }


    ### assign taxo
    taxo <- unique(taxo)

    if (length(taxo) == 0) { taxo <- 'Unassigned' }

    if (length(taxo) > 1) {

        stop('Taxonomy conflict in cluster ', x, ' : ', paste(taxo, collapse = ', '))

    }

    return(taxo)

})




# ------------------------------ simple termini ------------------------------ #
DF$simple_termini <- sapply(rownames(DF), function (x) {

    y <- t$genome_circ[ t$genome_id == x ]
    ifelse(y == '', 'noTR', 'TR')

})




# ------------------------------ simple source ------------------------------ #
DF$simple_source <- sapply(rownames(DF), function (x) {

    boo <- (t$source[ t$sp_repres == x ] == 'NL_vir')

    if (all(boo)) {

        S <- 'This study'

    } else if (any(boo)) {

        S <- 'Both'

    } else {

        S <- 'Databases'

    }

    return(S)

})




# ------------------------------ write table ------------------------------ #
write.table(
    DF,
    sep = '\t',
    quote = F,
    file = 'DB2_info.txt'
)




# ------------------------------ % NL_vir contigs represented by detected, TR, caudo vOTUs ------------------------------ #
sele <- rownames(DF)[ (DF$simple_termini == 'TR') & (DF$simple_taxo == 'Caudoviricetes') ]

n3 <- sum(t$source[ t$sp_repres %in% sele ] == 'NL_vir')

cat('\nNL_vir contigs belonging to detected, TR, caudo vOTUs:', n3, '(', round(n3/(n1+n2)*100, 2), '% )\n')




# ------------------------------ detected noTR vOTUs from databases ------------------------------ #
DF$full_taxo <- sapply(rownames(DF), function (x) {

    v <- t$taxo_refseq[ t$sp_repres == x ]
    v <- unique(v)
    if (length(v) > 1) { v <- v[ v != '' ] }
    paste(v, collapse = '|')

})

sele <- rownames(DF)[ (DF$simple_termini == 'noTR') & (DF$simple_source != 'This study') ]

write.table(
    DF[sele, ],
    sep = '\t',
    quote = F,
    file = 'DB2_noTR_from_databases.txt'
)




# ------------------------------ detected taxonomically assigned non-caudo vOTUs ------------------------------ #
sele <- rownames(DF)[ !(DF$simple_taxo %in% c('Unassigned', 'Caudoviricetes')) ]

write.table(
    DF[sele, ],
    sep = '\t',
    quote = F,
    file = 'DB2_taxo_non_caudo.txt'
)
