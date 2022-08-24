.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(rhmmer)
library(IRanges)
sessionInfo()




genome_id <- commandArgs(trailingOnly = T)[1]




# ------------------------------ input file ------------------------------ #
file1 <- paste0(genome_id, '_vs_Pfam.txt')

t <- read_domtblout(file1)

t <- as.data.frame(t, stringsAsFactors = F)




# ------------------------------ collect data per protein-profile pair ------------------------------ #
DF <- data.frame(NULL)

for (x in unique(t$domain_name)) {

    for (y in unique(t$query_accession)) {

        idx <- which((t$domain_name == x) & (t$query_accession == y))
        if (length(idx) == 0) { next }

        ir <- IRanges(start = t$env_from[idx], end = t$env_to[idx])
        cl <- reduce(ir)
        cl <- as.data.frame(cl)
        if (sum(cl$width) < 100) { next }

        df <- data.frame(
            protein_id = x,
            profile_id = y,
            profile_name = t$query_name[idx][1],
            protein_env_cov = sum(cl$width),
            protein_env_coo = paste(paste0(cl$start, '-', cl$end), collapse = ';'),
            stringsAsFactors = F
        )

        DF <- rbind(DF, df)

    }

}




# ------------------------------ for each protein select profile providing maximal coverage ------------------------------ #
sele <- c()

for (x in unique(DF$protein_id)) {

    idx <- which(DF$protein_id == x)
    sele <- c(sele, idx[ DF$protein_env_cov[idx] == max(DF$protein_env_cov[idx]) ][1])

}

DF <- DF[sele, ]




# ------------------------------ output file ------------------------------ #
file2 <- paste0(genome_id, '_annot.txt')

write.table(DF, sep = '\t', row.names = F, quote = F, file = file2)
