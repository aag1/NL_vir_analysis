.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(rhmmer)
sessionInfo()




# ------------------------------ viral refseq info ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_taxo.txt',
    sep = '\t',
    header = T,
    quote = '',
    fill = T,
    stringsAsFactors = F
)

L <- sapply(tab$genome_id, function (x) {

    v <- strsplit(tab[tab$genome_id == x, 'taxonomy'], ';')[[1]]
    if (length(v) == 0) { v <- '' }
    return(v)

})

t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/test_cenote/viral_refseq_aa_nt_ids.2.txt',
    sep = '\t',
    header = T,
    quote = '',
    fill = T,
    stringsAsFactors = F
)




# ------------------------------ markers list ------------------------------ #
Q <- list(
    Caudoviricetes = c('TerL_caudo', 'TerL_crAss', 'TerL_VOG00461'),
    Herpesviridae = 'Herpes_MCP_PF03122_full',
    Papillomaviridae = 'Papilloma_MCP_L1_VOG05075',
    Polyomaviridae = 'Polyoma_coat_PF00718_full',
    Adenoviridae = 'Adeno_Hexon_protein_VOG05391',
    Tectiliviricetes_without_Adenoviridae = c('Bam_Toil_MCP', 'FLiP_group_MCP', 'Gemmatimonas_MCP', 'PM2_MCP', 'PRD1_MCP', 'STIV_MCP', 'odin_MCP'),
    Nucleocytoviricota = c('GVOGm0003_MCP', 'GVOGm0022_RNAPS', 'GVOGm0023_RNAPL', 'GVOGm0054_PolB', 'GVOGm0172_TFIIB', 'GVOGm0461_TOPOII', 'GVOGm0760_VLTF3')
)




# ------------------------------ hmmer hits ------------------------------ #
DF <- read_domtblout('markers_vs_viral_refseq.txt')
DF <- as.data.frame(DF, stringsAsFactors = F)

DF$genome_id <- sapply(DF$domain_name, function (x) t$genome_id[ t$protein_id == x ])




# ------------------------------ NCLVD ViralRecall ------------------------------ #
d1 <- read.table(
    'nclvd_ViralRecall/nclvd_ViralRecall.summary.tsv',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ TerL motifs ------------------------------ #
d2 <- read.table(
    'viral_refseq_TerL_genomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ sensitivity & specificity ------------------------------ #
n_unassigned <- names(L)[ unlist(lapply(L, function (v) (length(v) == 1 & v[1] == '') | grepl('unclassified ((bacterial|archaeal) )*viruses', v[2]))) ]


H <- unique( DF$genome_id[ DF$query_name == Q[[ 'Herpesviridae' ]] ] )


for (g in names(Q)) {

    if (g == 'Tectiliviricetes_without_Adenoviridae') {

        n_all <- names(L)[unlist(lapply(L, function (v) ('Tectiliviricetes' %in% v) & (!('Adenoviridae' %in% v))))]

    } else {

        n_all <- names(L)[unlist(lapply(L, function (v) g %in% v))]

    }



    if (g == 'Caudoviricetes') {

        n_detected <- d2$genome_id

        n_detected <- n_detected[ !(n_detected %in% H) ]

    } else if (g == 'Nucleocytoviricota') {

        n_detected <- d1$replicon[ (d1$score >= 2) & (d1$markerhits != '-') ]

    } else {

        n_detected <- unique( DF$genome_id[ DF$query_name %in% Q[[g]] ] )

    }



    TP <- sum(n_detected %in% n_all)
    FN <- sum(!(n_all %in% n_detected))
    FP <- sum(!(n_detected %in% c(n_all, n_unassigned)))
    TN <- length(L) - length(unique(c(n_all, n_detected, n_unassigned)))

    sensitivity <- round(TP / (TP + FN) * 100, 2)
    specificity <- round(TN / (TN + FP) * 100, 2)

    cat('\n', g, ':\n')
    cat('Sensitivity:', sensitivity, '%\n')
    cat('Specificity:', specificity, '%\n\n')



    if (FP != 0) {

        sele <- n_detected[ !(n_detected %in% c(n_all, n_unassigned)) ]

        write.table(tab[tab$genome_id %in% sele, ], row.names = F, quote = F, sep = '\t', file = paste0(g, '_FP.txt'))

    }

    if (FN != 0) {

        sele <- n_all[ !(n_all %in% n_detected) ]

        write.table(tab[tab$genome_id %in% sele, ], row.names = F, quote = F, sep = '\t', file = paste0(g, '_FN.txt'))

    }

}
