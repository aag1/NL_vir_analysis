.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(seqinr)
sessionInfo()




option_list = list(make_option('--prefix'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)




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
L <- list()

for (n in c(11, 4, 15)) {

    # read gff
    f <- paste0(opt$prefix, '_c', n, '.gff')
    df <- read_gff(f, n)
    if (is.null(DF)) { DF <- df } else { DF <- merge(DF, df, by = 'genome_id', all = T) }


    # read fasta
    f <- paste0(opt$prefix, '_c', n, '.fasta')
    L[[paste0('c', n)]] <- read.fasta(f, seqtype = 'AA')

}

DF[ is.na(DF) ] <- 0




# ------------------------------ select code ------------------------------ #
DF$sele_code <- NA

for (i in 1:nrow(DF)) {

    # select code
    N <- DF$total_cscore_c11[i] * 1.1

    K <- 'c11'

    if (DF$total_cscore_c4[i] > N) { K <- 'c4' }

    if ((DF$total_cscore_c15[i] > N) & (DF$total_cscore_c15[i] > DF$total_cscore_c4[i])) { K <- 'c15' }

    DF$sele_code[i] <- K


    # extract proteome
    my_regexpr <- paste0('^', gsub('\\.', '\\\\.', DF$genome_id[i]), '_')

    idx <- grep(my_regexpr, names(L[[K]]))

    l <- L[[K]][idx]


    # proteome sequences
    write.fasta(
        sequences = l,
        names = names(l),
        nbchar = 80,
        file.out = paste0(opt$prefix, '_proteomes.fasta'),
        open = ifelse(i == 1, 'w', 'a')
    )


    # proteome info
    TAB <- data.frame(NULL)

    for (j in seq_along(l)) {

        A <- attributes(l[[j]])$Annot

        V <- strsplit(A, '#')[[1]]
        V[1] <- sub('^>', '', V[1])
        V <- sub('^ ', '', V)
        V <- sub(' $', '', V)

        tab <- data.frame(
            genome_id = sub('_[0-9]+_c[0-9]{1,2}$', '', V[1]),
            protein_id = V[1],
            start = as.numeric(V[2]),
            end = as.numeric(V[3]),
            strain = as.numeric(V[4]),
            stringsAsFactors = FALSE
        )
        TAB <- rbind(TAB, tab)

    }

    write.table(
        TAB,
        sep = '\t',
        quote = F,
        row.names = F,
        col.names = (i == 1),
        append = (i != 1),
        file = paste0(opt$prefix, '_proteomes.txt')
    )

}




# ------------------------------ write table ------------------------------ #
write.table(
    DF,
    sep = '\t',
    quote = F,
    row.names = F,
    file = paste0(opt$prefix, '_SeleCode.txt')
)
