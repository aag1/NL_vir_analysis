.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
sessionInfo()




# ------------------------------ VOG TerL profiles list ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/DATABASES/VOG_207/vog.annotations.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    comment.char = '',
    stringsAsFactors = FALSE
)

colnames(tab)[1] <- sub('^X\\.', '', colnames(tab)[1])

idx <- which(
    grepl('terminase', tab$ConsensusFunctionalDescription, ignore.case = TRUE) & 
    grepl('large', tab$ConsensusFunctionalDescription, ignore.case = TRUE)
)

tab <-  tab[idx, ]



d <- read.table(
    '/data/umcg-tifn/DATABASES/VOG_207/vog.lca.tsv',
    sep = '\t',
    quote = '',
    header = TRUE,
    comment.char = '',
    stringsAsFactors = FALSE
)

colnames(d)[1] <- sub('^X\\.', '', colnames(d)[1])



tab$LastCommonAncestor <- sapply(tab$GroupName, function (x) d$LastCommonAncestor[d$GroupName == x])



write.table(
    tab,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'vog_terl_profiles.txt'
)




# ------------------------------ VOG00461 taxonomic composition ------------------------------ #
V <- rownames(read.fasta('/data/umcg-tifn/DATABASES/VOG_207/vog_msa/VOG00461.msa')$ali)

V <- sub('^[0-9]+\\.', '', V)



d1 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_taxo.txt',
    sep = '\t',
    quote = '',
    fill = TRUE,
    header = TRUE,
    stringsAsFactors = FALSE
)

d2 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt',
    sep = '\t',
    quote = '',
    fill = TRUE,
    header = TRUE,
    stringsAsFactors = FALSE
)



DF <- data.frame(NULL)

for (protein_id in V) {

    genome_id <- d2$genome_id[d2$protein_id == protein_id]

    taxonomy <- d1$taxonomy[d1$genome_id == genome_id]

    df <- data.frame(protein_id, genome_id, taxonomy, stringsAsFactors = FALSE)
    DF <- rbind(DF, df)

}



write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'VOG00461_taxo.txt'
)
