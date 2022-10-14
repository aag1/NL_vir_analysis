.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(vioplot)
sessionInfo()




COL <- c(LLD = 'forestgreen', LLD2 = 'lightcoral', X300OB = 'royalblue1', IBD = 'goldenrod1', DEVoC_HA = 'grey75', DEVoC_all = 'grey75')

pdf('DB3_abundance_stability.pdf', height = 3.3, width = 7.5)

layout(matrix(1:2, ncol = 2), width = c(56, 44))




# ------------------------------ panel A ------------------------------ #
L <- list()
cat('\n')

for (x in names(COL)) {

    if (x %in% c('DEVoC_HA', 'DEVoC_all')) {

        f <- '../from_Peregrine/DEVoC_DB3_mapping_summary.txt'

    } else { f <- paste0('../from_Peregrine/', sub('^X', '', x), '_DB3_mapping_summary.txt') }


    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
    V <- setNames(t$pct_reads_mapped, t$sample_id)

    if (x == 'DEVoC_HA') {

       d <- read.table('../from_Peregrine/DEVoC_adult_cohort_PRJNA722819.tsv', sep = '\t', header = T, stringsAsFactors = F)
       V <- V[ d$run_accession[ grepl('^HA', d$sample_alias) ] ]

    }

    L[[x]] <- V


    cat(sub('^X', '', x), 'cohort, % sample reads mapped:\n')
    cat('min', round(min(V), 2), '%\n')
    cat('max', round(max(V), 2), '%\n')
    cat('mean', round(mean(V), 2), '%\n\n')

}


par(mar = c(2.5, 4.5, 1, 1), ps = 9)

vioplot(
    L,
    col = COL,
    names = NA,
    pchMed = 20,
    ylim = c(0, 90),
    frame.plot = F, xaxt = 'n', yaxt = 'n'
)

axis(side = 2, at = seq(0, 90, 30), las = 1)
mtext(side = 2, line = 2.5, text = 'Mapped reads, %', cex = 11/9)

N <- sub('^X', '', names(L))
N <- sub('^DEVoC_', 'DFV ', N)
N[N == 'LLD'] <- 'LLD\nbaseline'
N[N == 'LLD2'] <- 'LLD\nfollow-up'
text(x = seq_along(N) - 0.25, y = -12, labels = N, srt = 45, xpd = TRUE)

mtext('A', side = 3, line = -0.75, at = -1.25, cex = 24/9, xpd = TRUE)




# ------------------------------ panel B ------------------------------ #
load('DB3_permutation_test_data.RData')

par(mar = c(2.5, 5, 2.75, 0.5))

vioplot(
    V3, V1, V2,
    col = c('grey50', 'forestgreen', 'lightcoral'),
    names = NA,
    pchMed = 20,
    ylim = c(0.2, 1),
    las = 1, frame.plot = F, xaxt = 'n'
)

text(
    x = c(1, 2.15, 3.15),
    y = 0.1,
    labels = c('Intra-individual', 'LLD\nbaseline', 'LLD\nfollow-up'),
    xpd = TRUE
)

mtext(side = 2, line = 2.5, text = 'Bray-Curtis dissimilarity', cex = 11/9)

lines(x = c(1, 2), y = rep(1.05, 2), xpd = T)
text(x = 1.5, y = 1.08, labels = 'P-value < 0.0001', cex = 8/9, xpd = T)

lines(x = c(1, 3), y = rep(1.15, 2), xpd = T)
text(x = 2, y = 1.18, labels = 'P-value < 0.0001', cex = 8/9, xpd = T)

mtext('B', side = 3, line = 1, at = -0.75, cex = 24/9, xpd = TRUE)

dev.off()
