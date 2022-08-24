sessionInfo()




# ------------------------------ function ------------------------------ #
read_gff <- function (file, code) {

    t <- read.table(file, sep = '\t', header = F, stringsAsFactors = F)

    t$cscore <- as.numeric(lapply(strsplit(t$V9, ';'), function (v) sub('cscore=', '', v[9])))

    df <- data.frame(
        unique(t$V1),
        sapply(unique(t$V1), function (x) sum(t$cscore[t$V1 == x])),
        stringsAsFactors = F
    )

    colnames(df) <- c('genome_id', paste0('total_cscore_c', code))

    return(df)

}




# ------------------------------ read gff files ------------------------------ #
DF <- NULL

for (n in c(11, 4, 15)) {

    f <- paste0('CRASS_DB_cl_c', n, '.gff')

    df <- read_gff(f, n)

    if (is.null(DF)) { DF <- df } else { DF <- merge(DF, df, by = 'genome_id', all = T) }

}

DF[ is.na(DF) ] <- 0




# ------------------------------ select code ------------------------------ #
DF$sele_code <- 'c11'

for (i in 1:nrow(DF)) {

    N <- DF$total_cscore_c11[i] * 1.1

    if (DF$total_cscore_c4[i] > N) { DF$sele_code[i] <- 'c4' }

    if ((DF$total_cscore_c15[i] > N) & (DF$total_cscore_c15[i] > DF$total_cscore_c4[i])) { DF$sele_code[i] <- 'c15' }

}




# ------------------------------ compare with existing prediction ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)


DF$real_code <- tab[DF$genome_id, 'translation']

DF$real_code[ DF$real_code == 'TAG|q' ] <- 'c15'
DF$real_code[ DF$real_code == 'TGA|w' ] <- 'c4'

DF$success <- ifelse(DF$sele_code == DF$real_code, 'yes', 'no')


pct <- round(sum(DF$success == 'yes') / nrow(DF) * 100)
cat('\n\nCode prediction was successful in', sum(DF$success == 'yes'), 'cases out of', nrow(DF), '(', pct, '%)\n\n')


DF[DF$success == 'no', ]


write.table(DF, sep = '\t', row.names = F, quote = F, file = 'crAss378_predicted_code.txt')
