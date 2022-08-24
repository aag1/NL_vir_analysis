.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(RColorBrewer)
library(bio3d)
library(ggseqlogo)
library(ggplot2)
library(png)
sessionInfo()




# ------------------------------ data ------------------------------ #
ALI <- read.fasta('DB3_terL_msa_trimmed.fasta')$ali

M <- list(
    Walker_A = c(43, 51),
    Walker_B = c(128, 140),
    Nuclease_I = c(266, 275),
    Nuclease_II = c(321, 327),
    Nuclease_III = c(400, 410)
)

COL <- setNames(brewer.pal(5, 'Set3'), names(M))




# ------------------------------ conservation ------------------------------ #
conserv <- conserv(ALI, method = 'similarity', sub.matrix = 'blosum62')

# 11 columns window sliding with a 1 column step
w_center <- 6:(ncol(ALI) - 5)

w_conserv <- sapply(w_center, function (i) { sum(conserv[(i-5):(i+5)]) / 11 })




# ------------------------------ logo plots ------------------------------ #
logo_input <- function (coo, ali) {

    from <- coo[1]
    to <- coo[2]

    v <- c()
    for (s in rownames(ali)) {

        v <- c(v, paste(ali[s, from:to], collapse=''))

    }

    return(v)

}



for (x in names(M)) {

    V <- logo_input(M[[x]], ALI)

    png_file <- paste0(x, '.png')

    png(png_file, type = 'cairo')
    obj <- ggseqlogo(V) + coord_cartesian(ylim = c(0, 4.5)) 
    print(obj)
    dev.off()

}




# ------------------------------ final plot ------------------------------ #
pdf('DB3_terL_msa_trimmed.pdf', height = 4.5, width = 7.5)

layout(matrix(c(rep(1, 5), 2:6), nrow = 2, byrow = T), height = c(6, 4))



# conservation profile
par(mar = par()$mar + c(0, 1.5, 0, 0))

plot(
    NA,
    xlim = c(0, ncol(ALI)),
    ylim = range(w_conserv),
    xlab = 'MSA columns',
    ylab = 'Conservation',
    cex.lab = 1.25,
    xaxt = 'n', yaxt = 'n', bty = 'n'
)

for (x in names(M)) {

    rect(
        xleft = M[[x]][1],
        xright = M[[x]][2],
        ybottom = par()$usr[3],
        ytop = 0.2,
        col = COL[x],
        border = NA,
        xpd = T
    )

    text(
        gsub('_', ' ', x),
        x = sum(M[[x]]) / 2,
        y = 0.24,
        xpd = T
    )

}

lines(x = w_center, y = w_conserv)

axis(side = 1, at = seq(from = 0, to = 450, by = 50), xpd = T)

axis(side = 2, at = seq(from = -0.1, to = 0.2, by = 0.1), las = 1, xpd = T)



# motif logos
par(mar = rep(1, 4))

for (x in names(M)) {

    plot(
        NA,
        xlim = c(0, 1),
        ylim = c(0, 1),
        axes = F,
        ann  = F
    )

    png_file <- paste0(x, '.png')

	rasterImage(
        image   = readPNG(png_file),
        xleft   = -0.12,
        xright  = 1.05,
        ybottom = -0.28,
        ytop    = 1.1
    )

    box(col = COL[x], lwd = 2)

}

dev.off()
