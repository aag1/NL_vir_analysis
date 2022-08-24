plot_nt_content <- function (l, Xmax, Ymax, COL) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0, Ymax),
        xlab = '', ylab = 'Nucleotide content, %',
        xaxs = 'i', yaxs = 'i',
        xaxt = 'n',
        bty = 'n'
    )


    Xcoo <- 501 + 0:(length(l[[1]]) - 1) * 200


    for (n in c('A', 'T', 'G', 'C')) {

        lines(x = Xcoo, y = l[[n]], col = COL[n])

    }

}
