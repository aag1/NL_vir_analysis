sessionInfo()



tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = F
)

V <- read.table(
    'NL_vir_rrna_contigs.ids',
    header = F,
    stringsAsFactors = F
)[, 1]



idx <- which(tab$final_name %in% V)
tab[idx, ]



write.table(
    tab[-idx, ],
    sep = '\t',
    row.names = F,
    quote = F,
    file = 'NL_vir_genome_fragments.no_rrna.txt'
)
