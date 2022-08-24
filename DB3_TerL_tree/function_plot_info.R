sessionInfo()




plot_info <- function (DF, YLIM, tick_coef = 1, max_len = 360, pp_cex = 0.25) {

    # ------------------------------ genome length ------------------------------ #
    barplot(
        DF$repres_len,
        horiz = T,
        col = 'grey15', border = NA,
        xlim = c(0, max_len * 1000), ylim = YLIM - 0.5,
        xaxs = 'i',  yaxs = 'i',
        space = 0
    )

    text(
        'Length, kb',
        x = sum(par()$usr[1:2]) / 2,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )

    axis(side = 3, line = -1, at = c(0, max_len * 1000), labels = F)

    text(
        c(0, max_len),
        x = c(0, max_len * 1000),
        y = YLIM[2] * tick_coef,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ genome GC content ------------------------------ #
    barplot(
        DF$repres_pgc,
        horiz = T,
        col = 'grey15', border = NA,
        xlim = c(20, 70), ylim = YLIM - 0.5,
        xaxs = 'i',  yaxs = 'i',
        space = 0, xpd = F
    )

    text(
        'GC content, %',
        x = sum(par()$usr[1:2]) / 2,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )

    axis(side = 3, line = -1, at = c(20, 70), labels = F)

    text(
        c(20, 70),
        x = c(20, 70),
        y = YLIM[2] * tick_coef,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ terminal repeats ------------------------------ #
    plot(
        NA,
        xlim = c(0, 1), ylim = YLIM,
        xaxs = 'i',  yaxs = 'i',
        axes = F, ann = F
    )

    idx <- which(DF$repres_TR == 'ITR')
    points(x = rep(0.5, length(idx)), y = idx, pch = 16, cex = pp_cex)

    text(
        'ITR',
        x = sum(par()$usr[1:2]) / 2,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ genetic code ------------------------------ #
    plot(
        NA,
        xlim = c(0, 1), ylim = YLIM,
        xaxs = 'i',  yaxs = 'i',
        axes = F, ann = F
    )

    idx <- which(DF$repres_code == 'c15')
    points(x = rep(0.5, length(idx)), y = idx, pch = 16, cex = pp_cex, col = 'blue')

    idx <- which(DF$repres_code == 'c4')
    points(x = rep(0.5, length(idx)), y = idx, pch = 16, cex = pp_cex, col = 'red')

    text(
        'Alternative code',
        x = sum(par()$usr[1:2]) / 2,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ host phylum ------------------------------ #
    rle <- rle(DF$host_phylum_consensus)


    plot(
        NA,
        xlim = c(0, 1), ylim = YLIM,
        xaxs = 'i',  yaxs = 'i',
        axes = F, ann = F
    )


    for (i in seq_along(rle$values)) {

        if (rle$values[i] == '') { next }

        from <- ifelse(i==1, 1, sum(rle$lengths[1:(i-1)])+1)
        to <- from + rle$lengths[i] - 1

        rect(
            xleft = 0.35,
            xright = 0.65,
            ybottom = from - 0.5,
            ytop = to + 0.5,
            col = COL[ rle$values[i] ],
            border = NA
        )

    }


    text(
        'Host phylum',
        x = sum(par()$usr[1:2]) / 2,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ prophages in vOTU ------------------------------ #
    plot(
        NA,
        xlim = c(0, 1.25), ylim = YLIM,
        xaxs = 'i',  yaxs = 'i',
        axes = F, ann = F
    )

    idx <- which(DF$sp_prophage == 1)
    points(x = rep(0.5, length(idx)), y = idx, pch = 16, cex = pp_cex)

    text(
        'Prophage\n(Cenote-Taker2)',
        x = 0.5,
        y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
        srt = 90, adj = 0, xpd = T
    )




    # ------------------------------ plot % positive samples ------------------------------ #
    COL <- c(LLD = 'forestgreen', LLD2 = 'lightcoral', X300OB = 'royalblue1', IBD = 'goldenrod1', DEVoC_HA = 'grey75', DEVoC_all = 'grey75')

    for (k in names(COL)) {

        k_ <- sub('^X', '', k)
        x <- paste0('p_samples_', k_)

        barplot(
            DF[, x],
            horiz = T,
            col = COL[k], border = NA,
            xlim = c(0, 30), ylim = YLIM - 0.5,
            xaxs = 'i',  yaxs = 'i',
            space = 0
        )

        k_ <- sub('^DEVoC_', 'DFV ', k_)
        if (k_ == 'LLD2') { k_ <- 'LLD f/u' }
        text(
            k_,
            x = sum(par()$usr[1:2]) / 2,
            y = par()$usr[3] + diff(par()$usr[3:4]) * 1.035,
            srt = 90, adj = 0, xpd = T
        )

        axis(side = 3, line = -1, at = c(0, 30), labels = F)

        text(
            c(0, 30),
            x = c(0, 30),
            y = YLIM[2] * tick_coef,
            srt = 90, adj = 0, xpd = T
        )

    }

}
