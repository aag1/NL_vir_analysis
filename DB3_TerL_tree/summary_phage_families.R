.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(phangorn)
sessionInfo()




# -------------------------------------- read data -------------------------------------- #
ids3 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB3.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids3 <- rev(ids3)

ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids4 <- rev(ids4)

tab <- read.table('/data/umcg-tifn/NL_vir_analysis/cluster_genomes/DB1_info.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

tr <- read.tree('DB3_terL_tree.rooted.newick')




# -------------------------------------- build families matrix -------------------------------------- #
# refseq family
tab$family_refseq <- sapply(tab$taxo_refseq, function (x) {

    v <- strsplit(x, ';')[[1]]
    f <- grep('viridae$', v, value = T)
    if (length(f) == 0) { f <- '' } else { f <- f[1] }
    return(f)

})

tab$family_refseq[tab$family_refseq %in% c('Myoviridae', 'Podoviridae', 'Siphoviridae')] <- ''



# Benler_2020 family
tab$family_Benler_2020 <- sapply(tab$taxo_Benler_2020, function (x) {

    v <- strsplit(x, ',')[[1]]
    f <- grep('viridae$', v, value = T)
    if (length(f) == 0) { f <- '' } else { f <- f[1] }
    return(f)

})



# crAss-like family
tab$family_crAss <- ''
tab$family_crAss[tab$taxo_our_crAss != ''] <- 'crAss-like'



# list families
FAM <- unique(unlist(tab[, c('family_refseq', 'family_Benler_2020', 'family_crAss', 'taxo_BanfieldLab')]))
FAM <- FAM[FAM != '']



# create matrix
M <- matrix(
    0,
    nrow = length(ids3),
    ncol = length(FAM),
    dimnames = list(ids3, FAM)
)



# fill matrix
for (x in rownames(M)) {

    idx <- which(tab$sp_repres == x)

    fam <- unique(unlist(tab[idx, c('family_refseq', 'family_Benler_2020', 'family_crAss', 'taxo_BanfieldLab')]))
    fam <- fam[fam != '']

    M[x, fam] <- 1

}

sele <- apply(M, 2, function (v) any(v == 1))
M <- M[, sele]



# are there any DB3 vOTUs assigned to multiple families?
boo <- any(apply(M, 1, sum) > 1)
cat('\n\n\nAre there any DB3 vOTUs assigned to multiple families?', boo, '\n\n\n\n')



# write matrix
write.table(M, quote = F, sep = '\t', file = 'DB3_phage_families_matrix.txt')




# -------------------------------------- extend family assignments based on tree -------------------------------------- #
DF <- data.frame(
    phage_family = sapply(ids3, function (x) { idx <- which(M[x, ] != 0); ifelse(length(idx) == 0, '', colnames(M)[idx]) }),
    phage_family_ext = '',
    stringsAsFactors = F
)

rownames(DF) <- ids3



for (f in colnames(M)) {

    g <- rownames(M)[ M[, f] == 1 ]



    # family members --> MRCA --> MRCA descendants
    if (length(g) == 1) {

        mrca_desc <- g

    } else {

        mrca <- getMRCA(tr, tip = g)
        mrca_desc <- Descendants(tr, node = mrca, type = 'tips')[[1]]
        mrca_desc <- tr$tip.label[mrca_desc]

    }



    if (all(DF[mrca_desc, 'phage_family_ext'] == '')) {

        DF[mrca_desc, 'phage_family_ext'] <- f

    } else {

        stop('There are representatives of other families among the ', f, ' family MRCA descendants!')

    }

}



write.table(DF, quote = F, sep = '\t', file = 'DB3_phage_families_extended_via_mrca.txt')



cat('Initially,', sum(DF$phage_family != ''), 'vOTUs were taxonomically assigned:\n')
sort(table(DF$phage_family), decreasing = T)
cat('\n\n\n')



cat('After MRCA-based extension,', sum(DF$phage_family_ext != ''), 'vOTUs were taxonomically assigned:\n')
sort(table(DF$phage_family_ext), decreasing = T)
cat('\n\n\n')



DF <- DF[ids4, ]
boo <- apply(DF, 1, function (v) { (v[1] != '') | (v[2] != '') })
cat('DB4 phages assigned to a family:\n\n')
DF[boo, ]
cat('\n\n\n')
