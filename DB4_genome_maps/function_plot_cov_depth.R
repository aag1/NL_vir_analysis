plot_cov_depth <- function (l, Xmax, Ypower, COL) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(1, 10^Ypower),
        xlab = '', ylab = 'Mean read depth + 1',
        xaxs = 'i', yaxs = 'i',
        xaxt = 'n', yaxt = 'n', log = 'y',
        bty = 'n'
    )


    Xcoo <- 501 + 0:(length(l[[1]][[1]]) - 1) * 200


    for (k in names(l)) {

            invisible(lapply(l[[k]], function (v) {

                    lines(x = Xcoo, y = v + 1, col = adjustcolor(COL[k], alpha.f = 0.2), xpd = T)

            }))

    }


    axis(side = 2, at = 10^c(0:Ypower), labels = as.expression(sapply(0:Ypower, function (i) bquote(10^.(i)))))

}
