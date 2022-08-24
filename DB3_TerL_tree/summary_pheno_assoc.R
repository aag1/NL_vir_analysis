sessionInfo()



ids4 <- read.table('/data/umcg-tifn/NL_vir_analysis/DB3_DB4_info/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]
ids4 <- rev(ids4)

t1 <- read.table('/home/umcg-agulyaeva/NL_vir_analysis/DB3_TerL_tree/from_Gearshift/LLD_pheno_assoc.txt', sep = '\t', header = T, stringsAsFactors = F)

t2 <- read.table('/home/umcg-agulyaeva/NL_vir_analysis/DB3_TerL_tree/from_Gearshift/LLD_vs_IBD.txt', sep = '\t', header = T, stringsAsFactors = F)

t3 <- read.table('/home/umcg-agulyaeva/NL_vir_analysis/DB3_TerL_tree/from_Gearshift/LLD_vs_300OB.txt', sep = '\t', header = T, stringsAsFactors = F)



# pch 2  - overrepr. (w/o adj. for host)
# pch 24 - overrepr.
# pch 6  - underrepr. (w/o adj. for host)
# pch 25 - underrepr.



DF <- matrix(NA, nrow = length(ids4), ncol = 3, dimnames = list(ids4, c('lld_pheno', 'assoc_ibd', 'assoc_300ob')))
DF <- as.data.frame(DF, stringsAsFactors = F)



idx <- which(t1$FDR.glm1 < 0.05)
DF[t1$taxon[idx], 'lld_pheno'] <- t1$pheno[idx]



sele <- t2$taxon[ (t2$FDR.glm1 < 0.05) & (t2$LLD_positive_samples < t2$IBD_positive_samples) ]
DF[sele, 'assoc_ibd'] <- 2

sele <- t2$taxon[ (t2$FDR.glm2 < 0.05) & (t2$LLD_positive_samples < t2$IBD_positive_samples) ]
DF[sele, 'assoc_ibd'] <- 24

sele <- t2$taxon[ (t2$FDR.glm1 < 0.05) & (t2$LLD_positive_samples > t2$IBD_positive_samples) ]
DF[sele, 'assoc_ibd'] <- 6

sele <- t2$taxon[ (t2$FDR.glm2 < 0.05) & (t2$LLD_positive_samples > t2$IBD_positive_samples) ]
DF[sele, 'assoc_ibd'] <- 25



sele <- t3$taxon[ (t3$FDR.glm1 < 0.05) & (t3$LLD_positive_samples < t3$X300OB_positive_samples) ]
DF[sele, 'assoc_300ob'] <- 2

sele <- t3$taxon[ (t3$FDR.glm2 < 0.05) & (t3$LLD_positive_samples < t3$X300OB_positive_samples) ]
DF[sele, 'assoc_300ob'] <- 24

sele <- t3$taxon[ (t3$FDR.glm1 < 0.05) & (t3$LLD_positive_samples > t3$X300OB_positive_samples) ]
DF[sele, 'assoc_300ob'] <- 6

sele <- t3$taxon[ (t3$FDR.glm2 < 0.05) & (t3$LLD_positive_samples > t3$X300OB_positive_samples) ]
DF[sele, 'assoc_300ob'] <- 25



write.table(DF, quote = F, sep = '\t', file = 'DB4_pheno_assoc_pch.txt')
