# Analysis of *Caudoviricetes* phages with genome terminal repeats in four Dutch cohorts


## Folders

* __test_cenote, test_code_prediction, benchmark_markers__ - test procedures

* __run_cenote__ - run Cenote-Taker 2 for the 4 cohorts

* __rrna_NL_vir__ - detect rRNA genes in the virus-like contigs from the 4 cohorts

* __predict_proteomes__ - predict proteomes for the virus-like contigs from the 4 cohorts

* __taxo_NL_vir__ - assign taxonomy to the virus-like contigs from the 4 cohorts

* __get_databases_caudo__ - detect *Caudoviricetes* genomes with terminal repeats in 6 databases

* __cluster_genomes__ - dereplicate virus-like contigs from the 4 cohorts and genomes from the databases together (build DB0, DB1)

* __map_reads__ - map reads from the 4 cohorts to DB1

* __DEVoC_reads__ - map reads from the 254 Danish fecal viromes to DB1

* __DB2_info, DB3_DB4_info__ - build DB2, DB3, DB4

* __DB3_vConTACT2__ - DB3 gene-sharing with TerL-encoding viruses

* __DB3_TerL_tree__ - build & plot phylogenetic tree

* __DB3_abundance_stability__ - % sample reads mapped to DB3 & Bray-Curtis dissimilarity

* __DB3_tRNA_scan__ - tRNA detection in DB3 vOTU representatives

* __host_prediction__ - prophage-, CRISPR- and co-abundance-based host prediction

* __DB4_in_GenBank, DB4_in_BanfieldLab__ - look for phages similar to DB4 vOTU representatives

* __DB4_genome_info__ - nt content & GC skew for extended DB4

* __DB4_proteome_annot__ - functional annotation of proteins encoded by DB4 vOTU representatives

* __DB4_genome_maps__ - plot genome maps of DB4 vOTU representatives

* __DB4_pheno_assoc__ - associations between DB4 vOTUs prevalence and human phenotypes

* __workflow_diagram__ - plot analysis workflow

* __blastn_crAss_sensu_stricto__ - compare MGV-GENOME-0359371 to the first discovered crAssphage

* __DEVoC_LoVEphage__ - LoVEphage in the current study?


## Designations

* __DB0__ - 100,060 genomes & contigs pooled together

* __DB1__ - 30,461 vOTU representatives, result of DB0 dereplication

* __DB2__ - 15,196 detected vOTU representatives, result of read mapping to DB1

* __DB3__ - 1,899 DB2 vOTUs represented by *Caudoviricetes* genomes with terminal repeats

* __DB4__ - 54 DB4 vOTUs detected in > 5% samples in a Dutch cohort
