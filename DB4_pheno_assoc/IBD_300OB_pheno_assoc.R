.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(foreach)
sessionInfo()




# -------------------- read data --------------------
G <- read.table('../from_Peregrine/DB4.tree.ids', header = F, stringsAsFactors = F)[, 1]



### keys
key_lld <- read.table('../../LLD_LLD2_Info/LLD_GTMicrob.txt', sep = '\t', header = F, stringsAsFactors = F)
key_lld[, 2] <- paste0('X', key_lld[, 2])

key_ibd <- read.table('../../IBD_Info/rename_IBD.txt', sep = '\t', header = T, stringsAsFactors = F)

key_300ob <- read.table('../../300OB_Info/key_300OB.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)



### phages
phages <- list()

for (k in c('LLD', 'IBD', 'X300OB')) {

    file <- paste0('../from_Peregrine/', ifelse(k == 'X300OB', '300OB', k), '_DB1_abundance.txt')

    df <- read.table(file, row.names = 1, sep = '\t', header = T, stringsAsFactors = F)

    if (k == 'IBD') {
        z <- read.table('../../IBD_Info/selected_ibd_samples.txt', header = F, stringsAsFactors = F)[, 1]
        df <- df[, z]
    }

    df <- t(df)

    df <- df[, G]

    phages[[k]] <- df

}



### LLD phenotypes
pheno_lld <- data.frame(NULL)

pheno_files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = T)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(pheno_lld) == 0) { pheno_lld <- df } else { pheno_lld <- cbind(pheno_lld, df) }

}

rownames(pheno_lld) <- sapply(rownames(pheno_lld), function (x) key_lld[ key_lld[, 1] == x, 2 ])



### IBD phenotypes
pheno_ibd <- read.table(
    '../../IBD_Info/data_from_Renate/LLD_IBD_meta_201020.txt',
    quote = '',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_ibd) %in% key_ibd$Classic[key_ibd$Classic != 'XXXX'])
rownames(pheno_ibd)[idx] <- sapply(rownames(pheno_ibd)[idx], function (x) key_ibd[ key_ibd$Classic == x, 'old' ])
rownames(pheno_ibd)[rownames(pheno_ibd) == 'XXXX'] <- 'YYYY'

idx <- which(rownames(pheno_ibd) %in% key_lld[, 1])
rownames(pheno_ibd)[idx] <- sapply(rownames(pheno_ibd)[idx], function (x) key_lld[ key_lld[, 1] == x, 2 ])



### 300OB phenotypes
pheno_300ob_1 <- read.table(
    '../../300OB_Info/300OB_65phenotype.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_300ob_1) %in% key_300ob$ID)
rownames(pheno_300ob_1)[idx] <- sapply(rownames(pheno_300ob_1)[idx], function (x) key_300ob[ key_300ob$ID == x, 'G_id' ])

pheno_300ob_2 <- read.table(
    '../../300OB_Info/4tb8spry4b-1/300OB_metabolicSyndrome.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_300ob_2) %in% key_300ob$ID)
rownames(pheno_300ob_2)[idx] <- sapply(rownames(pheno_300ob_2)[idx], function (x) key_300ob[ key_300ob$ID == x, 'G_id' ])



### microbes
metaphlan <- read.table('../../crAss_analysis/from_Peregrine/LLD_LLD2_300OB_IBD_merged_abundance_table.txt', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

metaphlan <- metaphlan[, colnames(metaphlan) != 'NCBI_tax_id']

colnames(metaphlan) <- sub('_metaphlan$', '', colnames(metaphlan))

metaphlan <- t(metaphlan)

sc <- min(metaphlan[metaphlan > 0]) / 2    # small constant



### microbes-hosts
tab <- read.table('../from_Peregrine/DB4_phage_host_coAbundance_metaAnalysis_top_corr.txt', sep = '\t', header = T, stringsAsFactors = F)

hosts <- as.data.frame(log10(metaphlan[, unique(tab$host)] + sc))    # log transformation to get normally distributed data




# -------------------- LLD vs. IBD --------------------
DF1 <- rbind(phages[['LLD']], phages[['IBD']])

patient_cohort <- c(rep(0, nrow(phages[['LLD']])), rep(1, nrow(phages[['IBD']])))
names(patient_cohort) <- rownames(DF1)

DF2 <- cbind(data.frame(patient_cohort), pheno_ibd[rownames(DF1), c('AgeAtFecalSampling', 'Sex')])



OUT1 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = DF2, family = 'binomial'))

    data.frame(
        taxon = p,
        LLD_positive_samples = round(sum(phages[['LLD']][, p] > 0) / nrow(phages[['LLD']]) * 100, 2),
        IBD_positive_samples = round(sum(phages[['IBD']][, p] > 0) / nrow(phages[['IBD']]) * 100, 2),
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT2 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]
    h <- tab$host[ tab$phage == p ]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = cbind(DF2, hosts[rownames(DF1), h]), family = 'binomial'))

    data.frame(
        taxon = p,
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT <- data.frame(
    taxon = OUT1$taxon,

    LLD_positive_samples = OUT1$LLD_positive_samples,
    IBD_positive_samples = OUT1$IBD_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'BH'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'BH')
)

OUT <- OUT[order(OUT$FDR.glm1), ]



write.table(
    OUT,
    file = 'LLD_vs_IBD.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- IBD: CD vs. UC --------------------
Diagnosis <- setNames(pheno_ibd$DiagnosisCurrent, rownames(pheno_ibd))
Diagnosis <- Diagnosis[ rownames(phages[['IBD']]) ]
Diagnosis <- Diagnosis[ !is.na(Diagnosis) ]

ids_ibd_cd <- names(Diagnosis)[Diagnosis == 'CD']
ids_ibd_uc <- names(Diagnosis)[Diagnosis == 'UC']
IDS <- c(ids_ibd_cd, ids_ibd_uc)

Diagnosis <- Diagnosis[IDS]
Diagnosis <- ifelse(Diagnosis == 'CD', 0, 1)

DF1 <- phages[['IBD']][IDS, ]

DF2 <- cbind(
    data.frame(Diagnosis),
    pheno_ibd[IDS, c('AgeAtFecalSampling', 'Sex')]
)



OUT1 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = DF2, family = 'binomial'))

    data.frame(
        taxon = p,
        CD_positive_samples = round(sum(DF1[ids_ibd_cd, i] > 0) / length(ids_ibd_cd) * 100, 2),
        UC_positive_samples = round(sum(DF1[ids_ibd_uc, i] > 0) / length(ids_ibd_uc) * 100, 2),
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT2 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]
    h <- tab$host[ tab$phage == p ]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = cbind(DF2, hosts[rownames(DF1), h]), family = 'binomial'))

    data.frame(
        taxon = p,
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT <- data.frame(
    taxon = OUT1$taxon,

    CD_positive_samples = OUT1$CD_positive_samples,
    UC_positive_samples = OUT1$UC_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'BH'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'BH')
)

OUT <- OUT[order(OUT$FDR.glm1), ]



write.table(
    OUT,
    file = 'IBD_CD_vs_UC.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- IBD: exclusively colonic vs. ileal-inclusive disease location --------------------
Location <- setNames(pheno_ibd$DiseaseLocation, rownames(pheno_ibd))
Location <- Location[ rownames(phages[['IBD']]) ]
Location <- Location[ !is.na(Location) ]

ids_ibd_ec <- names(Location)[Location == 'colon']
ids_ibd_ii <- names(Location)[Location %in% c('ileum', 'both')]

IDS <- c(ids_ibd_ec, ids_ibd_ii)
Location <- Location[IDS]
Location <- ifelse(Location == 'colon', 0, 1)

DF1 <- phages[['IBD']][IDS, ]

DF2 <- cbind(
    data.frame(Location),
    pheno_ibd[IDS, c('AgeAtFecalSampling', 'Sex')]
)



OUT1 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = DF2, family = 'binomial'))

    data.frame(
        taxon = p,
        ec_positive_samples = round(sum(DF1[ids_ibd_ec, i] > 0) / length(ids_ibd_ec) * 100, 2),
        ii_positive_samples = round(sum(DF1[ids_ibd_ii, i] > 0) / length(ids_ibd_ii) * 100, 2),
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT2 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]
    h <- tab$host[ tab$phage == p ]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = cbind(DF2, hosts[rownames(DF1), h]), family = 'binomial'))

    data.frame(
        taxon = p,
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT <- data.frame(
    taxon = OUT1$taxon,

    ec_positive_samples = OUT1$ec_positive_samples,
    ii_positive_samples = OUT1$ii_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'BH'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'BH')
)

OUT <- OUT[order(OUT$FDR.glm1), ]



write.table(
    OUT,
    file = 'IBD_excl_colonic_vs_ileum_incl.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- LLD vs. 300OB --------------------
DF1 <- rbind(phages[['LLD']], phages[['X300OB']])

ids_lld   <- rownames(phages[['LLD']])
ids_300ob <- rownames(phages[['X300OB']])

DF2 <- data.frame(
    patient_cohort = c(
        rep(0, nrow(phages[['LLD']])),
        rep(1, nrow(phages[['X300OB']]))),
    Age = c(
        pheno_lld[ids_lld, 'antrop_age'],
        pheno_300ob_1[ids_300ob, 'age']),
    Gender = c(
        pheno_lld[ids_lld, 'antrop_gender.F1M2'],
        pheno_300ob_1[ids_300ob, 'sex'])
)

rownames(DF2) <- rownames(DF1)



OUT1 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = DF2, family = 'binomial'))

    data.frame(
        taxon = p,
        LLD_positive_samples = round(sum(phages[['LLD']][, p] > 0) / nrow(phages[['LLD']]) * 100, 2),
        X300OB_positive_samples = round(sum(phages[['X300OB']][, p] > 0) / nrow(phages[['X300OB']]) * 100, 2),
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT2 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]
    h <- tab$host[ tab$phage == p ]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = cbind(DF2, hosts[rownames(DF1), h]), family = 'binomial'))

    data.frame(
        taxon = p,
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT <- data.frame(
    taxon = OUT1$taxon,

    LLD_positive_samples = OUT1$LLD_positive_samples,
    X300OB_positive_samples = OUT1$X300OB_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'BH'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'BH')
)

OUT <- OUT[order(OUT$FDR.glm1), ]



write.table(
    OUT,
    file = 'LLD_vs_300OB.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- 300OB: with vs. without metabolic syndrome --------------------
DF1 <- phages[['X300OB']]

Syndrome <- setNames(pheno_300ob_2$MetS.NCEPsom, rownames(pheno_300ob_2))
Syndrome <- Syndrome[ ids_300ob ]
Syndrome <- ifelse(Syndrome >= 3, 1, 0)

ids_300ob_no <- names(Syndrome)[Syndrome == 0]
ids_300ob_sy <- names(Syndrome)[Syndrome == 1]

DF2 <- cbind(
    data.frame(Syndrome),
    pheno_300ob_1[ids_300ob, c('age', 'sex')]
)



OUT1 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = DF2, family = 'binomial'))

    data.frame(
        taxon = colnames(DF1)[i],
        no_positive_samples = round(sum(DF1[ids_300ob_no, i] > 0) / length(ids_300ob_no) * 100, 2),
        sy_positive_samples = round(sum(DF1[ids_300ob_sy, i] > 0) / length(ids_300ob_sy) * 100, 2),
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT2 <- foreach (i = 1:ncol(DF1), .combine = rbind) %do% {

    p <- colnames(DF1)[i]
    h <- tab$host[ tab$phage == p ]

    OBJ <- summary(glm(DF1[, i] > 0 ~ ., data = cbind(DF2, hosts[ids_300ob, h]), family = 'binomial'))

    data.frame(
        taxon = p,
        Beta = OBJ$coef[2,1],
        SE = OBJ$coef[2,2],
        Z = OBJ$coef[2,3],
        P = OBJ$coef[2,4]
    )

}

warnings()



OUT <- data.frame(
    taxon = OUT1$taxon,

    no_positive_samples = OUT1$no_positive_samples,
    sy_positive_samples = OUT1$sy_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'BH'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'BH')
)

OUT <- OUT[order(OUT$FDR.glm1), ]



write.table(
    OUT,
    file = '300OB_absence_vs_presence_metabolic_syndrome.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)
