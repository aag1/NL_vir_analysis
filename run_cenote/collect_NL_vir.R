.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



X <- c(
    'LLD' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD_scripts_logs_info/LLD_raw_reads_number_sele.txt',
    'LLD2' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD2_scripts_logs_info/LLD2_raw_reads_number.txt',
    '300OB' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/300OB_scripts_logs_info/300OB_raw_reads_number.txt',
    'IBD' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/IBD_scripts_logs_info/IBD_raw_reads_number_sele.txt'
)


t <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs.txt',
    header = T,
    sep = '\t',
    stringsAsFactors = FALSE
)



# ------------------------------ TABLE ------------------------------ #
TAB <- data.frame(NULL)

for (k in names(X)) {

    cat('Processing', k, 'cohort ...\n')

    samples <- read.table(X[k], sep = '\t', header = F, stringsAsFactors = F)[, 1]

    for (s in samples) {

        file <- paste0('final_combined_virus_sequences_', s, '_vir.fna')

        if (!file.exists(file)) { warning('File ', file, ' is missing!'); next() }

        L <- read.fasta(file, seqtype = 'DNA', forceDNAtolower = F)

        A <- unlist(lapply(L, function(x) attr(x, which = 'Annot')))

        DF <- data.frame(NULL)

        for (x in A) {

            x <- sub('^>', '', x)
            v <- strsplit(x, ' ')[[1]]
            if (length(v) == 3) { v <- v[c(1, 2, 2, 3)] }


            contig_len <- as.numeric(sub('^NODE_[0-9]+_length_([0-9]+)_cov_[0-9\\.]+$', '\\1', v[3]))


            if (v[2] != v[3]) {
                contig_fragm_coo <- v[2]
            } else {
                contig_fragm_coo <- paste0(0, '-', contig_len + 1)
            }


            idx <- which((t$cohort_name == k) & (t$sample_id == s) & (t$contig_id == v[3]))

            df <- data.frame(
                cohort = k,
                sample = s,
                contig = v[3],
                contig_fragm_coo,
                cenote_name = v[1],
                end_feature = ifelse(v[4] == 'no_end_feature', '', v[4]),
                crass_paper = ifelse(length(idx) == 0, '', t$contig_name[idx]),
                stringsAsFactors = F
            )

            DF <- rbind(DF, df)

        }

        DF$contig_fragm_len <- unlist(lapply(L, length))

        TAB <- rbind(TAB, DF)

    }

    sele <- which(TAB$cohort == k)
    TAB[sele, ] <- TAB[sample(sele), ]

}

warnings()

TAB$final_name <- sprintf('NL_vir%06d', c(1:nrow(TAB)))


write.table(
    TAB,
    sep = '\t',
    row.names = F,
    quote = F,
    file = 'NL_vir_genome_fragments.txt'
)



# ------------------------------ FASTA ------------------------------ #
for (s in unique(TAB$sample)) {

        file <- paste0('final_combined_virus_sequences_', s, '_vir.fna')

        L <- read.fasta(file, seqtype = 'DNA', forceDNAtolower = F)

        idx <- which(TAB$sample == s)

        write.fasta(
            sequences = L[ TAB$cenote_name[idx] ],
            names = TAB$final_name[idx],
            nbchar = 80,
            file.out = 'NL_vir_genome_fragments.fasta',
            open = ifelse(s == unique(TAB$sample)[1], 'w', 'a')
        )

}
