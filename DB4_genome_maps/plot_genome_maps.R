.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(TeachingDemos)
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_cov_depth.R')
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_genome_map.R')
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_nt_content.R')
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_nt_skew.R')
sessionInfo()




# ------------------------------ input data ------------------------------ #
ids <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', sep = '\t', header = F, stringsAsFactors = F)[, 1]



d <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)



### proteome coordinates
G <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)


for (db in c('MGV', 'GPD', 'GVD', 'DEVoC', 'HuVirDB', 'BanfieldLab')) {

    t <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_proteomes.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )

    G <- rbind(G, t)

}


t <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/Benler_2020_AA_proteomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

G <- rbind(G, t)



### proteome annotation
P <- list()

for (x in ids) {

    P[[x]] <- read.table(
        paste0('/data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot/', x, '_annot.txt'),
        sep = '\t',
        header = T,
        stringsAsFactors = F
    )

}



### tRNA coordinates
R <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/DB3_tRNA_scan/DB3_tRNAscan.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

R$genome_id <- sub(' $', '', R$genome_id)



### genome info
L <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_cov_depth.rds')

D <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_NT_CONTENT.rds')

GC_skew  <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_GC_SKEW.rds')
CGC_skew <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_CUMULATIVE_GC_SKEW.rds')

AT_skew  <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_AT_SKEW.rds')
CAT_skew <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_CUMULATIVE_AT_SKEW.rds')



### genome repeats
Q <- read.table('/data/umcg-tifn/NL_vir_analysis/DB4_proteome_annot/DB4_repeats.txt', sep = '\t', header = T, stringsAsFactors = F)



### colors
COL1 <- c(LLD = 'forestgreen', LLD2 = 'lightcoral', X300OB = 'royalblue1', IBD = 'goldenrod1')


COL2 <- c(A = 'darkred', T = 'darkseagreen', G = 'darkorchid', C = 'darksalmon')


v <- unique(unlist(lapply(P, function (t) t$profile_name)))

COL3 <- setNames(rep('grey70', length(v)), v)

x <- read.table(
    '/home/umcg-agulyaeva/NL_vir_analysis/DB4_proteome_annot/db4_detected_doms_function_manual.txt',
    sep = '\t',
    header = F,
    stringsAsFactors = F
)

sele <- x$V1[x$V2 %in% c('str', 'pack', 'assembly')]

COL3[names(COL3) %in% sele] <- 'steelblue1'

COL3[grepl('DNA_pol', names(COL3), ignore.case = T)] <- 'deeppink'

COL3[grepl('integr', names(COL3), ignore.case = T) | (names(COL3) == 'rve')] <- 'forestgreen'

COL3[names(COL3) == 'RVT_1'] <- 'darkorange2'




# ------------------------------ plot genome maps ------------------------------ #
pdf('DB4_genome_maps.pdf', height = 8.3, width = 11.7)

layout(
    matrix(1:5, nrow = 5),
    height = c(2, 3, 2, 1.5, 1.5)
)

par(
    las = 1,
    cex.lab = 1.25,
    mgp = par()$mgp + c(1, 0.5, 0.5),
    oma = c(4, 0, 1, 0),
    mar = par()$mar + c(-3, 2, -3, 5)
)

Xmax <- ceiling(max(d[ids, 'genome_len']) / 20000) * 20000



### for each genome ...
for (g in ids) {

    g_len <- d[g, 'genome_len']

    names(L[[g]])[ names(L[[g]]) == '300OB' ] <- 'X300OB'

    Xaxis <- seq(from = 0, to = ceiling(g_len / 20000) * 20000, by = 20000)



    ### coverage depth
    plot_cov_depth(L[[g]], Xmax, 3, COL1)

    axis(side = 1, at = Xaxis, labels = F)

    mtext(paste0(g, ' (', d[g, 'genome_circ'], ', code ', sub('^c', '', d[g, 'genome_code']), ')'), side = 3, line = 0)
    
    legend(
        x = max(Xaxis) + diff(par()$usr[1:2])*0.01,
        y = 10^2.5,
        legend = c('LLD', 'LLD f/u', '300OB', 'IBD'),
        col = sapply(COL1, function (x) adjustcolor(x, alpha.f = 0.2)),
        lwd = 1,
        xpd = T
    )



    ### genome map
    tab <- G[G$genome_id == g, 2:5]

    if (d[g, 'source'] == 'Benler_2020') {

        v <- strsplit(d[g, 'terL_coo'], ';')[[1]]

        strain <- ifelse(v[1] == 'f', 1, -1)
        start <- as.numeric(v[2])
        end <- as.numeric(v[3])

        te <- tab$protein_id[ (tab$start == start) & (tab$end == end) & (tab$strain == strain) ]

    } else { te <- d[g, 'terL_id'] }


    par(mar = par()$mar + c(0.5, 0, 7.5, 0))

    plot_genome_map(tab, te, P[[g]], R[R$genome_id == g, ], COL3, g_len, Xmax)


    q <- Q[Q$qseqid == g, ]

    if (nrow(q) > 0) {
        for (i in 1:nrow(q)) {

            y_coo <- 6.7 + 0.4 * (i - 1)
            l_type <- ifelse(q$sstrand[i] == 'plus', 1, 2)

            lines(x = c(q$from[i], q$to[i]), y = rep(y_coo + 0.2, 2), lty = l_type, col = 'darkorange2', xpd = T)
            lines(x = rep(q$from[i], 2), y = c(y_coo, y_coo + 0.2), lty = l_type, col = 'darkorange2', xpd = T)
            lines(x = rep(q$to[i], 2), y = c(y_coo, y_coo + 0.2), lty = l_type, col = 'darkorange2', xpd = T)

        }
    }



    ### nt content
    par(mar = par()$mar - c(0.5, 0, 7.5, 0))

    plot_nt_content(D[[g]], Xmax, 50, COL2)

    axis(side = 1, at = Xaxis, labels = F)

    legend(
        x = max(Xaxis) + diff(par()$usr[1:2])*0.01,
        y = 40,
        legend = names(COL2),
        col = COL2,
        lwd = 1,
        xpd = T
    )



    ### GC & AT skew
    r1 <- range(c(GC_skew[[g]], AT_skew[[g]]), na.rm = TRUE)

    r2 <- range(c(CGC_skew[[g]], CAT_skew[[g]]), na.rm = TRUE)


    plot_nt_skew(GC_skew[[g]], CGC_skew[[g]], 'GC', r1, r2, Xmax, max(Xaxis))

    axis(side = 1, at = Xaxis, labels = F)


    plot_nt_skew(AT_skew[[g]], CAT_skew[[g]], 'AT', r1, r2, Xmax, max(Xaxis))

    axis(side = 1, at = Xaxis, labels = Xaxis / 1000)

    mtext('Genome, kb', side = 1, line = par()$mgp[1], at = sum(range(Xaxis))/2, cex = 0.8, xpd = T)

}

dev.off()
