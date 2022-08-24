sessionInfo()



DF <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/DB2_info/DB2_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)



t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)



d <- data.frame(
    genome_id = rownames(DF)[ (DF$simple_termini == 'TR') & (DF$simple_taxo == 'Caudoviricetes') & (rownames(DF) != 'NL_vir005341') ],
    stringsAsFactors = F
)

d <- cbind(d, t[d$genome_id, c('terL_id', 'source')])



write.table(
    d,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB2_circ_caudo_terL.txt'
)
