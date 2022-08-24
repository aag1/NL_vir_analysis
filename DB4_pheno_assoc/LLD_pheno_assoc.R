.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(foreach)
sessionInfo()




# ------------------------------ read data ------------------------------ #
G <- read.table('../from_Peregrine/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]



key <- read.table('../../LLD_LLD2_Info/LLD_GTMicrob.txt', sep = '\t', header = F, stringsAsFactors = F)
key[, 2] <- paste0('X', key[, 2])



### phages
phages <- read.table('../from_Peregrine/LLD_DB1_abundance.txt', row.names = 1, sep = '\t', header = T, stringsAsFactors = F)

phages <- t(phages)

all(rownames(phages) %in% key[, 2])
rownames(phages) <- sapply(rownames(phages), function (x) key[ key[, 2] == x, 1 ], USE.NAMES = F)

phages <- phages[, G]



### phenotypes
pheno <- data.frame(NULL)

pheno_files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = T)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(pheno) == 0) { pheno <- df } else { pheno <- cbind(pheno, df) }

}

colnames(pheno)[colnames(pheno) == 'parasympathicolytic_inhaler'] <- 'parasympatholytic_inhaler'
colnames(pheno)[colnames(pheno) == 'antrop_heupomtrek'] <- 'antrop_hip.circ'
colnames(pheno)[colnames(pheno) == 'how_often_milk_or_sourmilk'] <- 'how_often_milk_or_buttermilk'
colnames(pheno)[colnames(pheno) == 'how_often_chocomilk_sweetend_milk_drinks'] <- 'how_often_chocomilk_sweetened_milk_drinks'

identical(sort(rownames(phages)), sort(rownames(pheno)))
pheno <- pheno[rownames(phages), ]



### microbes
metaphlan <- read.table('../../crAss_analysis/from_Peregrine/LLD_LLD2_300OB_IBD_merged_abundance_table.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

metaphlan <- metaphlan[, colnames(metaphlan) != 'NCBI_tax_id']

colnames(metaphlan) <- sub('_metaphlan$', '', colnames(metaphlan))

metaphlan <- t(metaphlan)

metaphlan <- metaphlan[(rownames(metaphlan) %in% key[, 2]), ]
rownames(metaphlan) <- sapply(rownames(metaphlan), function (x) key[ key[, 2] == x, 1 ], USE.NAMES = F)

sc <- min(metaphlan[metaphlan > 0]) / 2    # small constant

identical(sort(rownames(phages)), sort(rownames(metaphlan)))
metaphlan <- metaphlan[rownames(phages), ]



### microbes-hosts
tab <- read.table('../from_Peregrine/DB4_phage_host_coAbundance_metaAnalysis_top_corr.txt', sep = '\t', header = T, stringsAsFactors = F)

hosts <- as.data.frame(log10(metaphlan[, unique(tab$host)] + sc))   # log transformation to get normally distributed data




# ------------------------------ Logistic regression: phage presence/absence; adjusted for Age and Sex ------------------------------ #
result.glm1 <- foreach (i = 1:ncol(phages), .combine = rbind) %:%

  foreach (j = 1:ncol(pheno), .combine = rbind) %do% {

    if (colnames(pheno)[j] == 'antrop_age') {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    } else if (colnames(pheno)[j] == 'antrop_gender.F1M2') {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[3,1],
            SE = c1$coef[3,2],
            Z = c1$coef[3,3],
            P = c1$coef[3,4]
        )

    } else {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, j] + pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    }

  }

warnings()




# -------------------- Logistic regression: phage presence/absence; adjusted for Age, Sex and Host --------------------
result.glm2 <- foreach (i = 1:ncol(phages), .combine = rbind) %:%

  foreach (j = 1:ncol(pheno), .combine = rbind) %do% {

    p <- colnames(phages)[i]
    h <- tab$host[ tab$phage == p ]

    if (colnames(pheno)[j] == 'antrop_age') {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'] + hosts[, h], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    } else if (colnames(pheno)[j] == 'antrop_gender.F1M2') {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'] + hosts[, h], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[3,1],
            SE = c1$coef[3,2],
            Z = c1$coef[3,3],
            P = c1$coef[3,4]
        )

    } else {

        c1 <- summary(glm(phages[, i] > 0 ~ pheno[, j] + pheno[, 'antrop_age'] + pheno[, 'antrop_gender.F1M2'] + hosts[, h], family = 'binomial'))

        data.frame(
            pheno = colnames(pheno)[j],
            taxon = colnames(phages)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    }

  }

warnings()




# ------------------------------ multiple testing correction ------------------------------ #
result.glm1$FDR = p.adjust(result.glm1$P, method = 'BH')
result.glm2$FDR = p.adjust(result.glm2$P, method = 'BH')




# ------------------------------ combine results ------------------------------ #
combined = data.frame(
    pheno     = result.glm1$pheno,
    taxon     = result.glm1$taxon,

    Beta.glm1 = result.glm1$Beta,
    SE.glm1   = result.glm1$SE,
    Z.glm1    = result.glm1$Z,
    P.glm1    = result.glm1$P,
    FDR.glm1  = result.glm1$FDR,

    Beta.glm2 = result.glm2$Beta,
    SE.glm2   = result.glm2$SE,
    Z.glm2    = result.glm2$Z,
    P.glm2    = result.glm2$P,
    FDR.glm2  = result.glm2$FDR
)

combined = combined[order(combined$FDR.glm1), ]




# ------------------------------ write table ------------------------------ #
write.table(combined, sep = '\t', row.names = F, quote = F, file = 'LLD_pheno_assoc.txt')
