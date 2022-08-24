.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
sessionInfo()




t <- read.table(
    'DB4_cognate_GenBank_genera_level_with_pub_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

t <- t[(t$gb_id == 'KX501134.1') | grepl('^(Minot|Dzunkova)', t$authors), ]

t$frame_col <- ifelse(grepl('^Minot', t$authors), 'steelblue1', 'hotpink')
t$frame_col[t$gb_id == 'KX501134.1'] <- 'seagreen'




pdf('DB4_cognate_GenBank_dotplots.pdf', height = 11.7, width = 8.3)

layout(matrix(1:15, nrow = 5, byrow = T))

par(bty = 'n', mar = c(6, 1, 1, 1), oma = c(0, 0, 1.5, 0))


for (i in 1:nrow(t)) {

    cat('Working with the', t$nl_id[i], '&', t$gb_id[i], 'pair ...\n')


    ### matching words coo
    d <- read.table(
        paste0(t$nl_id[i], '___', t$gb_id[i], '.gff'),
        sep = '\t',
        header = F,
        stringsAsFactors = F
    )

    nl_ids <- which(d[, 1] == t$nl_id[i])
    gb_ids <- which(d[, 1] == t$gb_id[i])

    k <- data.frame(
        nl_from = d[nl_ids, 4],
        nl_to   = d[nl_ids, 5],
        gb_from = d[gb_ids, 4],
        gb_to   = d[gb_ids, 5],
        stringsAsFactors = F
    )


    ### plot canvas
    s <- abs(t$nl_len[i] - t$gb_len[i])

    if (t$nl_len[i] > t$gb_len[i]) {

        XLIM <- c(0 - s/2, t$gb_len[i] + s/2)
        YLIM <- c(0, t$nl_len[i])

    } else {

        XLIM <- c(0, t$gb_len[i])
        YLIM <- c(0 - s/2, t$nl_len[i] + s/2)

    }

    plot(
        NA,
        asp = 1,
        xlim = XLIM, ylim = YLIM,
        xaxs = 'i', yaxs = 'i',
        ann = F, axes = F, bty = 'n'
    )


    ### plot matching words
    apply(k, 1, function (v) lines(x = v[3:4], y = v[1:2], lwd = 0.1))


    ### plot annotation
    XLAB <- ifelse(
        t$gb_id[i] == 'KX501134.1',
        'contig71',
        grep('_', strsplit(t$gb_title[i], ' ')[[1]], value = T)
    )

    if (t$gb_rev[i] == 1) { XLAB <- paste(XLAB, 'RC') }

    YLAB <- sub('^(.+_NODE_[0-9]+)_length_[0-9]+_cov_[0-9\\.]+$', '\\1', t$nl_id[i])

    text(x = t$gb_len[i]/2, y = 0 - diff(par()$usr[3:4])*0.3, labels = XLAB, cex = 1.25, xpd = T)
    text(x = 0 - diff(par()$usr[1:2])*0.2, y = t$nl_len[i]/2, labels = YLAB, cex = 1.25, xpd = T, srt = 90)
    
    xtick <- floor(t$gb_len[i] / 1000)
    axis(side = 1, at = c(0, xtick*1000), labels = c(0, xtick), pos = 0, lwd = 0, lwd.ticks = 1)
    ytick <- floor(t$nl_len[i] / 1000)
    axis(side = 2, at = c(0, ytick*1000), labels = c(0, ytick), pos = 0, lwd = 0, lwd.ticks = 1)
    rect(xleft = 0, xright = t$gb_len[i], ybottom = 0, ytop = t$nl_len[i], border = t$frame_col[i], lwd = 1.5, xpd = T)

}


dev.off()
