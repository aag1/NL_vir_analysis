.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(vegan)
sessionInfo()




# ------------------------------ read data ------------------------------ #
M1 <- read.table(
    '../from_Peregrine/LLD_DB1_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = T,
    stringsAsFactors = F
)


M2 <- read.table(
    '../from_Peregrine/LLD2_DB1_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = T,
    stringsAsFactors = F
)


k <- read.table(
    '../../LLD_LLD2_Info/key_LLD_baseline_fup_338sample_pairs.txt',
    sep = '\t',
    row.names = 3,
    header = T,
    stringsAsFactors = F
)


db3 <- read.table(
    '../from_Peregrine/DB3.tree.ids',
    header = F,
    stringsAsFactors = F
)[, 1]

M1 <- M1[db3, ]
M2 <- M2[db3, ]




# ------------------------------ match LLD & LLD2 samples ------------------------------ #
MATCH <- sapply(colnames(M2), function (lld2_work_id) {

    lld_work_id <- k$barcode_baseline[ k$barcode_fup == paste0(lld2_work_id, '.bam') ]

    lld_work_id <- sub('^fece_', '', lld_work_id)

    lld_work_id <- paste0('X', lld_work_id)

    if (lld_work_id == 'ZZZZ') { lld_work_id <- 'QQQQ' }

    return(lld_work_id)

})

colnames(M2) <- paste0(MATCH, '_F')
M1 <- M1[, MATCH]




# ------------------------------ test stability of the overal phage composition ------------------------------ #
# Wilcoxon test on Bray-Curtis dissimilarity, with empiric P-value based on 10,000 permutations
# See Chen et al. 2021 (https://doi.org/10.1016/j.cell.2021.03.024)

m <- cbind(M1, M2)

m <- t(m)
m <- m[, colSums(m)>0]
m <- m[rowSums(m)>0, ]

bl_ids <- grep('_F', rownames(m), value = T, invert = T)    # baseline sample ids
fu_ids <- grep('_F', rownames(m), value = T, invert = F)    # follow-up sample ids

D <- vegdist(m, method = 'bray')
D <- as.matrix(D)


# 10,000 permutations
PERM_RES <- data.frame(NULL)

for (i in 1:10000) {

    d <- D
    n <- sample(rownames(D), nrow(D))
    rownames(d) <- n
    colnames(d) <- n


    d1 <- d[bl_ids, bl_ids]
    V1 <- d1[upper.tri(d1)]

    d2 <- d[fu_ids, fu_ids]
    V2 <- d2[upper.tri(d2)]

    sele <- which(paste0(bl_ids, '_F') %in% fu_ids)
    d3 <- d[bl_ids[sele], paste0(bl_ids[sele], '_F')]
    V3 <- diag(d3)


    res <- data.frame(
        pval_bl_fu = wilcox.test(V1, V2, alternative = 'two.sided', paired = F)$p.value,
        pval_in_bl = wilcox.test(V3, V1, alternative = 'two.sided', paired = F)$p.value,
        pval_in_fu = wilcox.test(V3, V2, alternative = 'two.sided', paired = F)$p.value,
        stringsAsFactors = F
    )
    PERM_RES <- rbind(PERM_RES, res)


    if (i %% 100 == 0) { cat('Iteration', i, '...\n') }

}


# real data
d1 <- D[bl_ids, bl_ids]
V1 <- d1[upper.tri(d1)]

d2 <- D[fu_ids, fu_ids]
V2 <- d2[upper.tri(d2)]

sele <- which(paste0(bl_ids, '_F') %in% fu_ids)
d3 <- D[bl_ids[sele], paste0(bl_ids[sele], '_F')]
V3 <- diag(d3)

REAL_RES <- data.frame(
    pval_bl_fu = wilcox.test(V1, V2, alternative = 'two.sided', paired = F)$p.value,
    pval_in_bl = wilcox.test(V3, V1, alternative = 'two.sided', paired = F)$p.value,
    pval_in_fu = wilcox.test(V3, V2, alternative = 'two.sided', paired = F)$p.value,
    stringsAsFactors = F
)


# empiric P-value
emp_pval_bl_fu <- sum( PERM_RES$pval_bl_fu < REAL_RES$pval_bl_fu ) / 10000
emp_pval_in_bl <- sum( PERM_RES$pval_in_bl < REAL_RES$pval_in_bl ) / 10000
emp_pval_in_fu <- sum( PERM_RES$pval_in_fu < REAL_RES$pval_in_fu ) / 10000


# save test results
save(V1, V2, V3, PERM_RES, REAL_RES, emp_pval_bl_fu, emp_pval_in_bl, emp_pval_in_fu, file = 'DB3_permutation_test_data.RData')
