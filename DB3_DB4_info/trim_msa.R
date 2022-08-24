.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
sessionInfo()




# -------------------------------------- read MSA -------------------------------------- #
msa <- read.fasta('TerL_msa.fasta')$ali




# -------------------------------------- DB2 terminal repeats, caudo, terL -------------------------------------- #
tab <- read.table('DB2_circ_caudo_terL.txt', sep = '\t', header = T, stringsAsFactors = F)

sele <- tab$terL_id[tab$terL_id != '']

msa <- msa[sele, ]

cat('\n\nvOTUs excluded due to undetected TerL in representative:', nrow(tab) - length(sele), '\n\n\n')




# -------------------------------------- key conserved TerL residues present -------------------------------------- #
# 10074, 2nd conserved acidic residue of Walker B motif
# 16944, 1st conserved acidic residue of nuclease motif I
# 17992, conserved acidic residue of nuclease motif II

IDX <- c()

for (i in c(10074, 16944, 17992)) {

    v <- msa[, i]
    idx <- which(!(v %in% c('D', 'E')))

    IDX <- c(IDX, idx)

    cat('Number of non-D/E in column', i, ':', length(idx), '(', round(length(idx)/nrow(msa)*100, 2), '%)\n')
    print(table(v[idx]))
    cat('\n\n')

}

IDX <- unique(IDX)
cat('vOTUs excluded due to non-D/E in key positions:', length(IDX), '\n\n')
msa <- msa[-IDX, ]

write.fasta(ids = rownames(msa), seqs = msa, file = 'DB3_terL_msa_full.fasta')




# -------------------------------------- MSA columns : less than 50% gaps -------------------------------------- #
pct_gap <- apply(msa, 2, function (v) sum(v == '-')) / nrow(msa) * 100

msa <- msa[, pct_gap < 50]




# -------------------------------------- write MSA -------------------------------------- #
write.fasta(ids = rownames(msa), seqs = msa, file = 'DB3_terL_msa_trimmed.fasta')

write.table(
    sapply(rownames(msa), function (x) tab$genome_id[ tab$terL_id == x ]),
    quote = F,
    row.names = F,
    col.names = F,
    file = 'DB3.ids'
)
