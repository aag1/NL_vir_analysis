library(igraph)
library(ggraph)
library(ggtext)



V <- read.table('proteomes.txt', sep = '\t', header = T, stringsAsFactors = F)[, 2:1]
V <- unique(V)
V$taxonomy[V$taxonomy == 'DB3'] <- 'CGTR1899'
V$taxonomy[V$genome %in% c('NC_055832.1', 'NC_055872.1', 'NC_002668.1', 'NC_008152.1', 'NC_047916.1')] <- 'CGTR1899'

E <- read.table('c1.ntw', header = F, stringsAsFactors = F)
colnames(E) <- c('from', 'to', 'weight')



g <- graph_from_data_frame(E, directed = F, vertices = V)



pdf('proteomes_vConTACT2_graph.pdf')

ggraph(graph = g, layout = 'stress') +

  geom_edge_link0(edge_colour = 'grey85', aes(edge_width = weight)) +

  geom_node_point(shape = 20, aes(color = taxonomy)) +

  scale_edge_width(range = c(0.1, 1), guide = 'none') +

  scale_color_manual(
    values = c('olivedrab', 'orange', 'steelblue1', 'firebrick1'),
    labels = c('CGTR1899', '*Caudoviricetes*', '*Herviviricetes*', '*Mirusviricota*'),
    name = '') +

  theme_graph(base_family = 'sans') +

  theme(
    legend.position = 'bottom',
    legend.background = element_rect(linetype = 'solid'),
    legend.text = element_markdown())

dev.off()
