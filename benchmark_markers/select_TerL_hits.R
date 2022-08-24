.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(rhmmer)
library(IRanges)
sessionInfo()




# ------------------------------ input data ------------------------------ #
option_list = list(
    make_option('--hmmerF'),
	make_option('--prefix'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


tab <- read_domtblout(opt$hmmerF)
tab <- as.data.frame(tab, stringsAsFactors = F)


tab <- tab[grepl('^TerL_', tab$query_name), ]


if (opt$prefix == 'viral_refseq') {

    d <- read.table(
        '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt',
        sep = '\t',
        header = T,
        quote = '',
        fill = T,
        stringsAsFactors = F
    )

    tab$genome_id <- sapply(tab$domain_name, function (x) d$genome_id[ d$protein_id == x ])

} else {

    tab$genome_id <- sub('_[0-9]+_c[0-9]+$', '', tab$domain_name)

}




# ------------------------------ hits table ------------------------------ #
MOTIFS <- list(

    TerL_caudo = c(
        WalkerA = 389,
        WalkerB = 483,
        NucleaseI = 612,
        NucleaseII = 667,
        NucleaseIII = 746
    ),

    TerL_crAss = c(
        WalkerA = 353,
        WalkerB = 465,
        NucleaseI = 681,
        NucleaseII = 744,
        NucleaseIII = 834
    ),

    TerL_VOG00461 = c(
        WalkerA = 340,
        WalkerB = 441,
        NucleaseI = 595,
        NucleaseII = 653,
        NucleaseIII = 732
    )

)


motif_names <- c('WalkerA', 'WalkerB', 'NucleaseI', 'NucleaseII', 'NucleaseIII')

for (m in motif_names) {

    tab[, m] <- 0

}


for (i in 1:nrow(tab)) {

    for (m in motif_names) {

        if ((tab$hmm_from[i] <= MOTIFS[[ tab$query_name[i] ]][m]) & (tab$hmm_to[i] >= MOTIFS[[ tab$query_name[i] ]][m])) {

            tab[i, m] <- 1

        }
    }
}


tab <- tab[order(tab$genome_id, tab$domain_name, tab$query_name), ]

tab <- tab[, c('genome_id', 'domain_name', 'query_name', 'sequence_evalue', motif_names, 'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'domain_len')]


write.table(
    tab,
    sep = '\t',
    quote = F,
    row.names = F,
    file = paste0(opt$prefix, '_TerL_hits.txt')
)




# ------------------------------ proteins table ------------------------------ #
TAB <- data.frame(NULL)


for (p in unique(tab$domain_name)) {

    idx <- which(tab$domain_name == p)

    protein_motifs <- apply(tab[idx, motif_names], 2, function (v) any(v == 1))
    protein_motifs <- ifelse(any(protein_motifs), paste(which(protein_motifs), collapse = ''), '')

    ir <- IRanges(start = tab$env_from[idx], end = tab$env_to[idx])
    cl <- reduce(ir)
    cl <- as.data.frame(cl)

    t <- data.frame(
        genome_id = tab$genome_id[idx][1],
        protein_id = p,
        protein_motifs,
        protein_cov = sum(cl$width),
        stringsAsFactors = F
    )
    TAB <- rbind(TAB, t)

}


write.table(
    TAB,
    sep = '\t',
    quote = F,
    row.names = F,
    file = paste0(opt$prefix, '_TerL_proteins.txt')
)




# ------------------------------ genomes table ------------------------------ #
DF <- data.frame(NULL)


for (g in unique(TAB$genome_id)) {

    idx <- which(TAB$genome_id == g)

    if (any(grepl('234', TAB$protein_motifs[idx]))) {

        idx <- idx[ grepl('234', TAB$protein_motifs[idx]) ]
        idx <- idx[ nchar(TAB$protein_motifs[idx]) == max(nchar(TAB$protein_motifs[idx])) ]
        idx <- idx[ TAB$protein_cov[idx] == max(TAB$protein_cov[idx]) ][1]

        df <- TAB[idx, c('genome_id', 'protein_id', 'protein_motifs')]
        DF <- rbind(DF, df)

    }
}


write.table(
    DF,
    sep = '\t',
    quote = F,
    row.names = F,
    file = paste0(opt$prefix, '_TerL_genomes.txt')
)
