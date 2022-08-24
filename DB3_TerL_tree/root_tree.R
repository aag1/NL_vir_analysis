.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(phangorn)
sessionInfo()



# ------------------------------ read data ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

tr <- read.tree('DB3_terL_tree.treefile')



# ------------------------------ rename tips, root  ------------------------------ #
tr$tip.label <- sapply(tr$tip.label, function (x) tab$genome_id[ tab$terL_id == x ])

tr <- midpoint(tr)

tr <- rotate(tr, node = 1900)

write.tree(tr, file = 'DB3_terL_tree.rooted.newick')



# ------------------------------ basic plot ------------------------------ #
tr <- read.tree('DB3_terL_tree.rooted.newick')

pdf('DB3_terL_tree.rooted.pdf')

plot(tr, show.tip.label = F, x.lim = 8)

tiplabels(tr$tip.label, cex = 0.25, frame = 'none', col = 'seagreen', adj = 0)

nodelabels(cex = 0.25, frame = 'none', col = 'seagreen', adj = 1)

tiplabels(tr$tip.label, cex = 0.25, frame = 'none', col = 'seagreen', adj = 0)

add.scale.bar()

dev.off()
