sessionInfo()



tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.txt', 
    header = TRUE, 
    sep = '\t', 
    quote = '', 
    fill = TRUE, 
    stringsAsFactors = FALSE
)



idx <- which(tab$genome_id %in% tab$protein_id)   # polyprotein instead of genome

for (i in idx) {

    j <- which(tab$protein_id == tab$genome_id[i])

    tab$genome_id[i] <- tab$genome_id[j]

}



write.table(
    tab,
    quote = F,
    sep = '\t',
    row.names = F,
    file = '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt'
)
