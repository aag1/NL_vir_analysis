V <- read.table('proteomes.txt', sep = '\t', header = T, stringsAsFactors = F)[, 2:1]
V <- unique(V)
V$taxonomy[V$taxonomy == 'DB3'] <- 'CGTR1899'
V$taxonomy[V$genome %in% c('NC_055832.1', 'NC_055872.1', 'NC_002668.1', 'NC_008152.1', 'NC_047916.1')] <- 'CGTR1899'

E <- read.table('c1.ntw', header = F, stringsAsFactors = F)
colnames(E) <- c('from', 'to', 'weight')



E$from_taxo <- sapply(E$from, function (x) V$taxonomy[V$genome == x])
E$to_taxo <- sapply(E$to, function (x) V$taxonomy[V$genome == x])




sum((E$from_taxo == 'CGTR1899') & !(E$to_taxo %in% c('CGTR1899', 'Caudoviricetes')))    # 0
sum((E$to_taxo == 'CGTR1899') & !(E$from_taxo %in% c('CGTR1899', 'Caudoviricetes')))    # 0



x1 <- unique(E$from[ (E$from_taxo == 'CGTR1899') & (E$to_taxo == 'Caudoviricetes') ])
x2 <- unique(E$to[ (E$to_taxo == 'CGTR1899') & (E$from_taxo == 'Caudoviricetes') ])
length(unique(x1, x2)) / sum(V$taxonomy == 'CGTR1899') * 100    # 51.8694



for (x in unique(V$taxonomy)) { cat(x, '\t', sum(V$taxonomy == x), '\n') }
# Mirusviricota 	 111 
# Herviviricetes 	 92 
# Caudoviricetes 	 4162 
# CGTR1899 	 1899 
