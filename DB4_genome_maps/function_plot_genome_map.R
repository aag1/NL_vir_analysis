# library TeachingDemos must be loaded!


plot_genome_map <- function (tab, terL_id, dom, trna, col, len, Xmax) {

    ### canvas
    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(6.5, 0.5),
        xaxs = 'i', yaxs = 'i', 
        axes = FALSE, ann = FALSE, bty = 'n'
    )



 	### ORFs
    DOM_LABELS <- c()

	for (i in 1:nrow(tab)) {

		from   <- tab[i, 2]
		to     <- tab[i, 3]
        strand <- tab[i, 4]


        if (strand == 1) {
            frame <- ((from - 1) %% 3) + 1
            Y <- frame
        }
        if (strand == -1) {
            frame <- ((from - 1) %% 3) + 1
            Y <- frame + 3
        }


        # grey ORF rectangle
		rect(
            xleft = from,
            xright = to,
            ybottom = Y + 0.5,
            ytop = Y - 0.5,
            col = 'grey98',
            border = NA
        )


        # protein domains
        idx <- which(dom$protein_id == tab[i, 1])

        if (length(idx) == 1) {

            d <- dom$profile_name[idx]
            s <- strsplit(dom$protein_env_coo[idx], ';')[[1]]
            l <- lapply(strsplit(s, '-'), as.numeric)
            c1 <- unlist(lapply(l, function (v) v[1]))
            c2 <- unlist(lapply(l, function (v) v[2]))

            if (strand == 1) {
                d_from <- (from - 1) + c1*3 - 2
                d_to <- (from - 1) + c2*3
            }

            if (strand == -1) {
                d_from <- (to + 1) - c2*3
                d_to <- (to + 1) - c1*3 + 2
            }

            rect(
                xleft = d_from,
                xright = d_to,
                ybottom = Y + 0.5,
                ytop = Y - 0.5,
                col = col[d],
                border = NA
            )

            DOM_LABELS <- c(DOM_LABELS, setNames((min(d_from) + max(d_to)) / 2, d))

        }


        # ORF border
		rect(
            xleft = from,
            xright = to,
            ybottom = Y + 0.5,
            ytop = Y - 0.5,
            col = NA,
            border = 'black'
        )

    }



    ### tRNA
    if (nrow(trna) > 0) {

	   for (i in 1:nrow(trna)) {

           from   <- min(trna$tRNA_begin[i], trna$tRNA_end[i])
           to     <- max(trna$tRNA_begin[i], trna$tRNA_end[i])
           strand <- ifelse(trna$tRNA_begin[i] < trna$tRNA_end[i], 1, -1)


           if (strand == 1) {
               frame <- ((from - 1) %% 3) + 1
               Y <- frame
           }
           if (strand == -1) {
               frame <- ((from - 1) %% 3) + 1
               Y <- frame + 3
           }


           rect(
               xleft = from,
               xright = to,
               ybottom = Y + 0.5,
               ytop = Y - 0.5,
               col = 'firebrick',
               border = 'firebrick'
           )

           name <- paste0('tRNA_', trna$tRNA_type[i], '_', trna$anticodon[i])
           DOM_LABELS <- c(DOM_LABELS, setNames((from + to) / 2, name))

       }

    }

    V <- grep('^tRNA_', unique(names(DOM_LABELS)), value = T)
    col <- c(col, setNames(rep('firebrick', length(V)), V))



	### genome
	rect(
        xleft = 0,
        xright = len, 
        ybottom = 6.5,
        ytop = 0.5,
        col = NA,
        border = 'black'
    )

	text(
        x = -0.01 * Xmax,
        y = 1:6,
        labels = c('+1', '+2', '+3', '-1', '-2', '-3'),
        adj = 1,
        xpd = TRUE
    )



    ### highlight TerL
    i <- which(tab[, 1] == terL_id)

    if (tab[i, 4] == 1) { Y <- ((tab[i, 2] - 1) %% 3) + 1 } else { Y <- ((tab[i, 2] - 1) %% 3) + 4 }

    rect(
        xleft = tab[i, 2],
        xright = tab[i, 3],
        ybottom = Y + 0.5,
        ytop = Y - 0.5,
        col = NA,
        border = 'red'
    )



    ### domain labels
    K <- spread.labs(DOM_LABELS, mindiff = Xmax * 0.01, min = 0, max = Xmax, maxiter = 10^8)

    text(
        names(DOM_LABELS),
        x = K,
        y = rep(-1, length(K)),
        col = col[ names(DOM_LABELS) ],
        cex = 0.8,
        srt = 90,
        adj = 0,
        xpd = T
    )

    for (i in seq_along(K)) {

        lines(
            x = c(DOM_LABELS[i], K[i]),
            y = c(0.1, -0.6),
            col = col[ names(DOM_LABELS)[i] ],
            lwd = 0.5,
            xpd = T
        )

    }

}
