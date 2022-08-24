sessionInfo()



TAB <- data.frame(NULL)

for (i in 0:9) {

    tab <- read.table(
        paste0('DB2_circ_', i, '-scores.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )

    if (nrow(TAB) == 0) { TAB <- tab } else { TAB <- rbind(TAB, tab) }

}



write.table(
    TAB,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'DB2_circ_PlasX_score.txt'
)
