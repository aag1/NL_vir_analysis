sessionInfo()



outF <- commandArgs(trailingOnly = T)[1]



DF <- data.frame(NULL)

for (f in list.files(pattern = '_tRNA_genes.txt$')) {

    if (file.info(f)$size == 0) { next }

    df <- read.table(f, sep = '\t', skip = 3, header = F, stringsAsFactors = F)
    df[, 2] <- NULL

    colnames(df) <- c('genome_id', 'tRNA_begin', 'tRNA_end', 'tRNA_type', 'anticodon', 'intron_begin', 'intron_end', 'score', 'note')

    DF <- rbind(DF, df)

}



write.table(DF, sep = '\t', row.names = F, quote = F, file = outF)
