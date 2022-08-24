try(source('/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/link_phage_host.R'))
sessionInfo()




# ------------------------------ Spacer database Shmakov et al. 2017 ------------------------------ #
# read data
tab <- read.table(
    'blastn_spacers_Shmakov2017_vs_DB3.spacer80match.out',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)

V <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'nident', 'qlen')
colnames(tab) <- V


# q_organism
tab$q_organism <- unlist(lapply(strsplit(tab$qaccver, '_'), function (v) v[1]))


# phage-host links
DF1 <- LINK(tab)




# ------------------------------ CRISPR-Cas++ spacer database 20210121 ------------------------------ #
# read data
t <- read.table(
    'blastn_spacers_CRISPRCasdb_vs_DB3.spacer80match.out',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(t) <- V

d1 <- read.table(
    '/data/umcg-tifn/DATABASES/CRISPRCasdb/20210121_spacer_34_safe_names.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


# q_organism
d1 <- d1[d1$spacer_id %in% unique(t$qaccver), ]

d2 <- data.frame(NULL)

for (i in 1:nrow(d1)) {

    x <- data.frame(spacer_id = d1$spacer_id[i], q_organism = strsplit(d1$genome_id[i], '\\|')[[1]], stringsAsFactors = FALSE)
    d2 <- rbind(d2, x)

}

tab <- merge(
        x = t,
        y = d2,
        by.x = 'qaccver',
        by.y = 'spacer_id',
        sort = FALSE
)


# phage-host links
DF2 <- LINK(tab)




# ------------------------------ output table ------------------------------ #
TAB <- rbind(DF1, DF2)

TAB <- aggregate(TAB$spacer, by = list(TAB$phage, TAB$host), function (v) paste(v, collapse = ';'))
colnames(TAB) <- c('phage', 'host', 'spacer')

write.table(
    TAB,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'DB3_crispr_links.txt'
)
