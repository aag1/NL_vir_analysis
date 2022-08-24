.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(rhmmer)
sessionInfo()




# ------------------------------ read data ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/rrna_NL_vir/NL_vir_genome_fragments.no_rrna.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

t <- read.table(
    'NL_vir_genome_fragments_circ.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ markers list ------------------------------ #
Q <- list(
    Caudoviricetes = c('TerL_caudo', 'TerL_crAss', 'TerL_VOG00461'),
    Herpesviridae = 'Herpes_MCP_PF03122_full',
    Adenoviridae = 'Adeno_Hexon_protein_VOG05391',
    Papillomaviridae = 'Papilloma_MCP_L1_VOG05075',
    Polyomaviridae = 'Polyoma_coat_PF00718_full',
    Tectiliviricetes = c('Bam_Toil_MCP', 'FLiP_group_MCP', 'Gemmatimonas_MCP', 'PM2_MCP', 'PRD1_MCP', 'STIV_MCP', 'odin_MCP'),
    Nucleocytoviricota = c('GVOGm0003_MCP', 'GVOGm0022_RNAPS', 'GVOGm0023_RNAPL', 'GVOGm0054_PolB', 'GVOGm0172_TFIIB', 'GVOGm0461_TOPOII', 'GVOGm0760_VLTF3')
)




# ------------------------------ hmmer hits ------------------------------ #
DF <- read_domtblout('markers_vs_NL_vir.txt')
DF <- as.data.frame(DF, stringsAsFactors = F)

DF$genome_id <- sub('_[0-9]+_c[0-9]+$', '', DF$domain_name)




# ------------------------------ NCLVD ViralRecall ------------------------------ #
d1 <- read.table(
    'nclvd_ViralRecall/nclvd_ViralRecall.summary.tsv',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ TerL motifs ------------------------------ #
d2 <- read.table(
    'NL_vir_TerL_genomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ terminal repeats ------------------------------ #
tab$contig_fragm_circ <- ''

for (x in unique(t$seq_id)) {

    idx1 <- which(tab$final_name == x)

    idx2 <- which(t$seq_id == x)

    tab$contig_fragm_circ[idx1] <- paste(t$term_rep_type[idx2], collapse = ';')


    # Note sequences identical to reverse complement, they may have emerged as follows:
    # 5'-ATGC-3'   --->   5'-ATGC-3'  +  5'-GCAT-3'   --->   5'-ATGCGCAT-3'
    # 3'-TACG-5'                                             3'-TACGCGTA-5'
    idx3 <- which((t$seq_id == x) & (t$term_rep_type == 'ITR'))

    if (length(idx3) == 1) {

        if (tab$contig_fragm_len[idx1] == t$term_rep_len[idx3]) {

            cat(tab$final_name[idx3], 'is identical to its reverse complement!\n')

        }

    }

}




# ------------------------------ taxonomic assignment ------------------------------ #
tab$contig_taxo <- ''


H <- unique( DF$genome_id[ DF$query_name == Q[[ 'Herpesviridae' ]] ] )


for (g in names(Q)) {

    # determine taxo
    if (g == 'Caudoviricetes') {

        detected <- d2$genome_id

        detected <- detected[ !(detected %in% H) ]

    } else if (g == 'Nucleocytoviricota') {

        detected <- d1$replicon[ (d1$score >= 2) & (d1$markerhits != '-') ]

    } else {

        detected <- unique( DF$genome_id[ DF$query_name %in% Q[[g]] ] )

    }


    # check for conflicts & assign taxo
    if (length(detected) > 0) {

        idx <- which(tab$final_name %in% detected)

        if (any(tab$contig_taxo[idx] != '')) {

            b <- data.frame(
                contid_name = tab$final_name[idx],
                taxo_1 = tab$contig_taxo[idx],
                taxo_2 = g,
                stringsAsFactors = F
            )

            cat('\nTaxo assignent conflict!\n'); print(b); cat('\n\n')

        }

        tab$contig_taxo[idx] <- g

    }

}


write.table(
    tab,
    sep = '\t',
    quote = F,
    row.names = F,
    file = 'NL_vir_genome_fragments.no_rrna.taxo.txt'
)
