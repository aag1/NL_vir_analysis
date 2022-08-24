.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()




sp_id <- commandArgs(trailingOnly = TRUE)[1]




# ------------------------------ read data ------------------------------ #
t1 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/run_cenote/NL_vir_genome_fragments.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

t2 <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)




# ------------------------------ get host fragments ------------------------------ #
SEQ <- list()

sele <- t2$genome_id[ (t2$sp_repres == sp_id) & (t2$source == 'NL_vir') & (t2$prophage_cenote == 1) ]

for (nm in sele) {

    idx <- which(t1$final_name == nm)


    chrt <- t1$cohort[idx]
    smpl <- t1$sample[idx]
    ctg  <- t1$contig[idx]
    len  <- as.numeric(sub('^NODE_[0-9]+_length_([0-9]+)_cov_[0-9\\.]+$', '\\1', ctg))
    coo  <- as.numeric(strsplit(t1$contig_fragm_coo[idx], '-')[[1]])


    file <- paste0('/data/umcg-tifn/', chrt, '_assembly/', smpl, '/', smpl, '_contigs.fasta')
    L <- read.fasta(file, seqtype = 'DNA', forceDNAtolower = FALSE)


    if (coo[1] != 0) {

        from <- 1
        to <- coo[1]
        if ((to - from + 1) >= 1000) {
            SEQ[[ paste0(ctg, '_', from, '_', to) ]] <- L[[ ctg ]][ from:to ]
        }

    }


    if (coo[2] != (len + 1)) {

        from <- coo[2]
        to <- len
        if ((to - from + 1) >= 1000) {
            SEQ[[ paste0(ctg, '_', from, '_', to) ]] <- L[[ ctg ]][ from:to ]
        }

    }

}




# ------------------------------ write data ------------------------------ #
if (length(SEQ) > 100) {

    SEQ <- SEQ[ sample(seq_along(SEQ), size = 100) ]

}

write.fasta(
    sequences = SEQ,
    names = names(SEQ),
    nbchar = 80,
    file.out = paste0(sp_id, '_host_fragm.fasta')
)
