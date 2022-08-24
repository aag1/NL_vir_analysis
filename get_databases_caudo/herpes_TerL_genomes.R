.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(rhmmer)
sessionInfo()



for (db in c('MGV', 'GPD', 'GVD', 'DEVoC', 'HuVirDB', 'BanfieldLab')) {

    file1 <- paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/markers_vs_', db, '_5kb_cMAG.txt')
    DF <- read_domtblout(file1)
    DF <- as.data.frame(DF, stringsAsFactors = F)
    DF$genome_id <- sub('_[0-9]+_c[0-9]+$', '', DF$domain_name)


    file2 <- paste0('/data/umcg-tifn/NL_vir_analysis/get_databases_caudo/', db, '_5kb_cMAG_TerL_genomes.txt')
    tab <- read.table(file2, sep = '\t', header = T, stringsAsFactors = F)


    H <- unique( DF$genome_id[ DF$query_name == 'Herpes_MCP_PF03122_full' ] )
    print(any(tab$genome_id %in% H))

}
