.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(RColorBrewer)
source('/home/umcg-agulyaeva/NL_vir_analysis/DB3_TerL_tree/function_plot_info.R')
sessionInfo()




# ------------------------------ read data ------------------------------ #
tr <- read.tree('DB3_terL_tree.rooted.newick')


DF <- read.table('DB3_data4plot.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)


M <- read.table('DB3_phage_families_matrix.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)
colnames(M)[colnames(M) == 'crAss.like'] <- 'crAss-like'


H <- read.table('DB3_predicted_host_phyla.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)
DF$host_phylum_consensus <- H[rownames(DF), 'host_phylum_consensus']


P <- read.table('DB4_pheno_assoc_pch.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)


v <- H$host_phylum_consensus
n <- names(sort(table(v[v != '']), decreasing = T))
COL <- setNames(brewer.pal(length(n), 'Paired'), n)


ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]


A <- read.table('DB4_pheno_assoc_pch.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)


t <- read.table('/data/umcg-tifn/NL_vir_analysis/DB4_in_GenBank/DB4_cognate_GenBank_genera_level_with_pub_info.txt', sep = '\t', header = T, stringsAsFactors = F)




# ------------------------------ DB3: plot tree ------------------------------ #
pdf('DB3_TerL_tree.pdf', height = 8.75, width = 7.5)

layout(matrix(1:14, ncol = 14), width = c(1.75, 1.35, rep(0.48, 5), 0.6, rep(0.5, 6)))

par(oma = c(0, 0, 0, 1), mar = c(0, 1, 9.5, 0))

plot(tr, show.tip.label = F)

add.scale.bar(length = 1)

YLIM = par()$usr[3:4]

legend(
    'topleft', inset = c(0, -0.15), xpd = T,
    title = 'Code',
    legend = c('15', '4'),
    col = c('blue', 'red'),
    pch = 16,
    ncol = 1,
    cex = 0.75
)

legend(
    'topleft', inset = c(0.34, -0.15), xpd = T,
    title = 'Host phylum',
    legend = names(COL),
    fill = COL,
    ncol = 1,
    cex = 0.75
)




# ------------------------------ DB3: plot taxonomy matrix ------------------------------ #
par(mar = c(0, 0.6, 9.5, 0.6))

plot(
    NA,
    xlim = c(0, ncol(M)+1),
    ylim = YLIM, yaxs = 'i',
    axes = F, ann = F
)

for (i in 1:ncol(M)) {

    lines(x = rep(i, 2), y = c(-10, Ntip(tr) + 10), col = 'grey95')

    idx <- which(M[, i] == 1)
    f <- colnames(M)[i]

    points(
        x = rep(i, length(idx)),
        y = idx,
        cex = 0.25,
        pch = 16
    )

    text(
        labels = f,
        x = i,
        y = par()$usr[4],
        font = ifelse(f %in% c('Autographiviridae', 'Drexlerviridae'), 3, 1),
        srt = 90,
        adj = 0,
        xpd = T
    )

}




# ------------------------------ DB3: plot info ------------------------------ #
plot_info(DF, YLIM)

mtext('Positive samples, %', side = 3, line = 7, adj = 1.5, cex = 0.8)

dev.off()




# ------------------------------ DB4: plot tree ------------------------------ #
tr <- keep.tip(tr, tip = ids4)

DF <- DF[tr$tip.label, ]

COL <- COL[ names(COL) %in% unique(DF$host_phylum_consensus) ]



tip_color <- setNames(rep('black', Ntip(tr)), tr$tip.label)
tip_color[c('uvig_191910', 'MGV-GENOME-0359371')] <- 'firebrick' # crAss-like
tip_color['MGV-GENOME-0351725'] <- 'seagreen'                    # Flandersviridae
tip_color['MGV-GENOME-0305083'] <- 'blueviolet'                  # Faecalibacterium phage FP_Lagaffe

L <- list()
L[[1]] <- which(tr$tip.label %in% t$nl_id[ grepl('^Minot', t$authors) ])
L[[2]] <- which(tr$tip.label %in% t$nl_id[ grepl('^Dzunkova', t$authors) ])
col <- c('steelblue1', 'hotpink')

tr$tip.label <- sapply(tr$tip.label, function (x) sub('^(.+_NODE_[0-9]+)_length_[0-9]+_cov_[0-9\\.]+$', '\\1', x))



pdf('DB4_TerL_tree.pdf', height = 8.75, width = 7.5)

layout(matrix(1:14, ncol = 14), width = c(2.5, 1, rep(0.49, 2), rep(0.38, 3), 0.48, rep(0.5, 6)))

par(oma = c(0, 0, 0, 1), mar = c(0, 1, 9.5, 0))



plot(tr, x.lim = 20, show.tip.label = F)

add.scale.bar(length = 1)

YLIM = par()$usr[3:4]



tiplabels(
    text = tr$tip.label,
    tip = 1:Ntip(tr),
    col = tip_color,
    frame = 'n',
    cex = 0.75,
    adj = 0, offset = 0.2,
    xpd = T
)



for (n in seq_along(L)) {

    for (i in L[[n]]) {

        tiplabels(
            tip = i,
            col = col[n],
            offset = 0.5 + 0.75 * strwidth(tr$tip.label[i]),
            pch = 16,
            cex = 1.25,
            adj = 1
        )

    }

}



legend(
    'topleft', inset = c(0, -0.15), xpd = T,
    title = 'Code',
    legend = '15',
    col = 'blue',
    pch = 16
)

legend(
    'topleft', inset = c(0.34, -0.15), xpd = T,
    title = 'Host phylum',
    legend = names(COL),
    fill = COL,
    ncol = 1
)




# ------------------------------ DB4: plot associations ------------------------------ #
idx <- which(!is.na(A$lld_pheno))
A$lld_pheno[A$lld_pheno == 'antrop_age'] <- 'age'
A$lld_pheno[A$lld_pheno == 'BioMK_ChromograninA_log'] <- 'CgA'
text(
    A$lld_pheno[idx],
    x = par()$usr[2],
    y = idx,
    font = 2,
    cex = 0.75,
    adj = 1, xpd = T
)


par(mar = c(0, 0.6, 9.5, 0.6))

plot(
    NA,
    xlim = c(0.5, 2.5), ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    axes = F, ann = F
)


idx <- which(!is.na(A$assoc_300ob))
pch <- A$assoc_300ob[idx]
points(
    x = rep(1, length(idx)),
    y = idx + ifelse(pch %in% c(6, 25), 0.1, 0),
    pch = pch,
    col = 'royalblue1',
    bg = 'royalblue1'
)


idx <- which(!is.na(A$assoc_ibd))
pch <- A$assoc_ibd[idx]
points(
    x = rep(2, length(idx)),
    y = idx + ifelse(pch %in% c(6, 25), 0.1, 0),
    pch = pch,
    col = 'goldenrod1',
    bg = 'goldenrod1'
)


text(
    'Differential\nabundance',
    x = 1.5,
    y = par()$usr[4] + diff(par()$usr[3:4])*0.124,
    cex = 1.2,
    xpd = T
)

text(
    c('300OB', 'IBD'),
    x = c(1, 2),
    y = par()$usr[4],
    srt = 90, adj = 0, xpd = T
)




# ------------------------------ DB4: plot info ------------------------------ #
plot_info(DF, YLIM, tick_coef = 0.99, max_len = 150, pp_cex = 0.8)

mtext('Positive samples, %', side = 3, line = 7, adj = 1.5, cex = 0.8)

dev.off()
