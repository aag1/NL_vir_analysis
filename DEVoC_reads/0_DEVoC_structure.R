sessionInfo()



# https://doi.org/10.1128/mSystems.00382-21
# MicrobLiver project
# Chronic liver disease adults (CLD?): 52 individuals
# Healthy adults (HA?): 52 individuals
# Healthy minors (HM?): 50 individuals
# Obese minors (OM?): 50 individuals, 2 samples per individual (baseline and 1-year follow-up)



# Pediatric cohort: https://www.ebi.ac.uk/ena/browser/view/PRJNA723467
# Adult cohort: https://www.ebi.ac.uk/ena/browser/view/PRJNA722819
# 'Show Column Selection' --> add 'sample_alias', Download report: TSV



for (k in c('pediatric', 'adult')) {

    if (k == 'pediatric') { f <- 'DEVoC_pediatric_cohort_PRJNA723467.tsv' }
    if (k == 'adult') { f <- 'DEVoC_adult_cohort_PRJNA722819.tsv' }
    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
    if (nrow(t) > length(unique(t$sample_alias))) { stop('Multiple runs per sample!') }


    cat('\n\n')
    cat(k, 'cohort:', nrow(t), 'samples\n')


    for (x in c('CLD', 'HA', 'HM', 'OM')) {

        r <- paste0('^', x)
        n <- sum(grepl(r, t$sample_alias))
        cat(k, 'cohort:', n, x, 'samples\n')

    }

}
