.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(rhmmer)
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_cov_depth.R')
source('/home/umcg-agulyaeva/crAss_analysis/genome_maps/function_plot_genome_map.R')
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_nt_content.R')
source('/home/umcg-agulyaeva/NL_vir_analysis/DB4_genome_maps/function_plot_nt_skew.R')
sessionInfo()




# ------------------------------ input data ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    row.names = 1,
    stringsAsFactors = F
)



G <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/predict_proteomes/NL_vir_genome_fragments_proteomes.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

H <- read_domtblout('/data/umcg-tifn/NL_vir_analysis/taxo_NL_vir/markers_vs_NL_vir.txt')
H <- as.data.frame(H, stringsAsFactors = F)



L <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_cov_depth.rds')

D <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_NT_CONTENT.rds')

GC_skew  <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_GC_SKEW.rds')
CGC_skew <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_CUMULATIVE_GC_SKEW.rds')

AT_skew  <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_AT_SKEW.rds')
CAT_skew <- readRDS('/data/umcg-tifn/NL_vir_analysis/DB4_genome_info/DB4_CUMULATIVE_AT_SKEW.rds')



COL1 <- c(LLD = 'forestgreen', LLD2 = 'lightcoral', X300OB = 'royalblue1', IBD = 'goldenrod1')

COL2 <- c(A = 'darkred', T = 'darkseagreen', G = 'darkorchid', C = 'darksalmon')




# ------------------------------ plot genome map ------------------------------ #
pdf('NL_vir005341_genome_map.pdf', height = 8.3, width = 11.7)

layout(matrix(1:5, nrow = 5))

par(
    las = 1,
    cex.lab = 1.25,
    mgp = par()$mgp + c(1, 0.5, 0.5),
    oma = c(4, 0, 1, 0),
    mar = par()$mar + c(-3, 2, -3, 5)
)



p <- 'NL_vir005341'

p_len <- tab[p, 'genome_len']

names(L[[p]])[ names(L[[p]]) == '300OB' ] <- 'X300OB'

Xmax <- ceiling(p_len / 20000) * 20000

Xaxis <- seq(from = 0, to = Xmax, by = 20000)



### coverage depth
plot_cov_depth(L[[p]], Xmax, 3, COL1)

axis(side = 1, at = Xaxis, labels = F)

mtext(paste0(p, ' (', tab[p, 'genome_circ'], ')'), side = 3, line = 0)
    
legend(
    x = max(Xaxis) + diff(par()$usr[1:2])*0.01,
    y = 10^2.5,
    legend = c('LLD', 'LLD f/u', '300OB', 'IBD'),
    col = sapply(COL1, function (x) adjustcolor(x, alpha.f = 0.2)),
    lwd = 1,
    xpd = T
)



### genome map
g <- G[G$genome_id == p, 2:5]

terl <- tab[p, 'terL_id']
idx <- which(H$domain_name == terl)
coo <- paste0(H$env_from[idx], '-', H$env_to[idx])
coo <- as.data.frame(coo)
rownames(coo) <- terl
colnames(coo) <- 'TerL'

plot_genome_map(g, coo, setNames('red', 'TerL'), p_len, Xmax)



### nt content
plot_nt_content(D[[p]], Xmax, 50, COL2)

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
r1 <- range(c(GC_skew[[p]], AT_skew[[p]]), na.rm = TRUE)

r2 <- range(c(CGC_skew[[p]], CAT_skew[[p]]), na.rm = TRUE)


plot_nt_skew(GC_skew[[p]], CGC_skew[[p]], 'GC', r1, r2, Xmax, max(Xaxis))

axis(side = 1, at = Xaxis, labels = F)


plot_nt_skew(AT_skew[[p]], CAT_skew[[p]], 'AT', r1, r2, Xmax, max(Xaxis))

axis(side = 1, at = Xaxis, labels = Xaxis / 1000)

mtext('Genome, kb', side = 1, line = par()$mgp[1], at = sum(range(Xaxis))/2, cex = 0.8, xpd = T)

dev.off()
