plot_nt_skew <- function (skew, c_skew, pair, r1, r2, Xmax, c_axis_pos) {

    Xcoo <- 501 + 0:(length(skew) - 1) * 200


    plot(
        x = Xcoo,
        y = skew,
        type = 'l',
        col = 'blue',
        xlim = c(0, Xmax),
        ylim = c(floor(r1[1]/0.2), ceiling(r1[2]/0.2))*0.2,
        xaxs = 'i', yaxs = 'i',
        axes = FALSE, ann = FALSE, bty = 'n'
    )

    N <- diff(par()$usr[1:2])

    axis(side = 2, line = par()$mgp[3], col = 'blue', col.axis = 'blue')

    mtext(paste0(pair, ' skew'), side = 2, line = par()$mgp[1], cex = 0.8, las = 0, col = 'blue')


    par(new = TRUE)


    plot(
        x = Xcoo,
        y = c_skew,
        type = 'l',
        col = 'red',
        xlim = c(0, Xmax),
        ylim = c(floor(r2[1]/20), ceiling(r2[2]/20))*20,
        xaxs = 'i', yaxs = 'i',
        axes = FALSE, ann = FALSE, bty = 'n'
    )

    axis(side = 4, line = par()$mgp[3], col = 'red', col.axis = 'red', pos = c_axis_pos + N*0.01)

    text(paste0('C', pair, ' skew'), x = c_axis_pos + N*0.06, y = sum(par()$usr[3:4])/2, srt = 90, cex = 1.25, xpd = T, col = 'red')

}
