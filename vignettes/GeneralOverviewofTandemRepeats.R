## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>", fig.width=6)
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(123)

## ----sessioninfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE------
base_path <- "/home/matteo/polybox/MSc_ACLS/master_thesis/TRALResultAnalysis/"

## ----Housekeeping, message=FALSE, warning=FALSE, collapse=TRUE-----------
library(TRALResultAnalysis)
tr_crcfavorable_path <- paste0(base_path, "inst/extdata/TRs_favorable_proteins_CRC_sp.tsv")
tr_crcunfavorable_path <- paste0(base_path, "inst/extdata/TRs_unfavorable_proteins_CRC_sp.tsv")
tr_wnt_path <- paste0(base_path, "inst/extdata/TRs_Wnt_proteins_CRC_sp.tsv")
dest_file_sp <- paste0(base_path, "data/swissprot_human.tsv")
dest_file_kin <- paste0(base_path, "data/swissprot_human_kinome.tsv")

## ----load data sets, eval=TRUE, echo=FALSE, message=TRUE, warning=TRUE----
tr_fav <- load_tr_annotations(tr_crcfavorable_path)
tr_unfav <- load_tr_annotations(tr_crcunfavorable_path)
tr_wnt <- load_tr_annotations(tr_wnt_path)

## ----Download sp_all if it doesnt exist----------------------------------
sp_url <- "https://www.uniprot.org/uniprot/?query=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20reviewed:yes&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,virus%20hosts,encodedon,database(Pfam),interactor,comment(ABSORPTION),feature(ACTIVE%20SITE),comment(ACTIVITY%20REGULATION),feature(BINDING%20SITE),feature(CALCIUM%20BIND),comment(CATALYTIC%20ACTIVITY),comment(COFACTOR),feature(DNA%20BINDING),ec,comment(FUNCTION),comment(KINETICS),feature(METAL%20BINDING),feature(NP%20BIND),comment(PATHWAY),comment(PH%20DEPENDENCE),comment(REDOX%20POTENTIAL),rhea-id,feature(SITE),comment(TEMPERATURE%20DEPENDENCE)&sort=organism"

# Download the swissprot file only if it doesn't already exist.
# Uncomment this, if you want to use the most recent available data!
if(!file.exists(dest_file_sp)){
    download.file(sp_url, destfile = dest_file_sp)
}

sp_all_fav <- load_swissprot(dest_file_sp, tr_fav)
sp_all_unfav <- load_swissprot(dest_file_sp, tr_unfav)
sp_all_wnt <- load_swissprot(dest_file_sp, tr_wnt)

tr_unfav_sp <- merge(x = tr_unfav, y = sp_all_unfav, by = "ID", all.x = TRUE)
tr_fav_sp <- merge(x = tr_fav, y = sp_all_fav, by = "ID", all.x = TRUE)
tr_wnt_sp <- merge(x = tr_wnt, y = sp_all_wnt, by = "ID", all.x = TRUE)

## ----fraction of TR containing proteins----------------------------------
# CRC favorable
length(unique(tr_fav_sp$ID))/403 # through unique() we ensure to count those proteins with >1 TR only once.

# CRC unfavorable
length(unique(tr_unfav_sp$ID))/286 

# Wnt pathway
length(unique(tr_wnt_sp$ID))
length(unique(tr_wnt_sp$ID))/644 

## ----warning=FALSE-------------------------------------------------------
TR_location(
  rbind(tr_fav_sp, tr_unfav_sp),
  byTRtype = TRUE)

## ----warning=FALSE-------------------------------------------------------
TR_location(tr_fav_sp, byTRtype = TRUE)
TR_location(tr_unfav_sp, byTRtype = TRUE)

## ----warning=FALSE-------------------------------------------------------
TR_location(
   rbind(tr_fav_sp, tr_unfav_sp, tr_wnt_sp),
   byTRtype = TRUE)

## ------------------------------------------------------------------------
# combine all TR from the three groups
tr_all <- rbind(tr_fav, tr_unfav, tr_wnt)

AAfreq_all <- AAfreq_in_TR(tr_all) 
AAfreq_all[base::order(AAfreq_all$aa_ratio, decreasing = TRUE),]

AAfreq_all_homo <- AAfreq_in_TR(tr_all[which(tr_all$l_type == "homo"),]) 
AAfreq_all_homo[base::order(AAfreq_all_homo$aa_ratio, decreasing = TRUE),]

AAfreq_all[which(AAfreq_all$aa == "L"),]
AAfreq_all_homo[which(AAfreq_all_homo$aa == "L"),]

## ----TODO is this necessary----------------------------------------------
# (AAfreq_fav <- AAfreq_in_TR(tr_fav))
# (AAfreq_unfav <- AAfreq_in_TR(tr_unfav))
# cor.test(AAfreq_fav$aa_ratio, AAfreq_unfav$aa_ratio, method = "spearman")

## ------------------------------------------------------------------------
AAratio_vs_Disorderpropensity(tr_all, plot_title = "CRC & Wnt-Pathway Associated Proteins")

## ------------------------------------------------------------------------
AAratio_vs_Disorderpropensity(sp_overall = TRUE)

## ----fav in wnt-pathway--------------------------------------------------
# CRC favorable proteins in Wnt pathway
tr_fav_sp$protein_name[which(tr_fav_sp$ID %in% tr_wnt_sp$ID)]

## ----unfav in wnt-pathway------------------------------------------------
# CRC unfavorable proteins in Wnt pathway
tr_unfav_sp$protein_name[which(tr_unfav_sp$ID %in% tr_wnt_sp$ID)]

## ----many TRs------------------------------------------------------------
# CRC favorable proteins
table(table(tr_fav$ID))

# CRC unfavorable proteins
table(table(tr_unfav_sp$ID))

# Wnt-pathway proteins
table(table(tr_wnt_sp$ID))

## ----multiple TR regions in unfavorable prots I--------------------------
protein_id_by_number_of_TR(tr_unfav_sp, 11)
eleven_TRs <- tr_unfav_sp[which(tr_unfav_sp$ID == protein_id_by_number_of_TR(tr_unfav_sp, 11)[[1]]),]
unique(eleven_TRs$protein_name)
unique(eleven_TRs$prot_function)

## ----multiple TR regions in unfavorable prots II-------------------------
protein_id_by_number_of_TR(tr_unfav_sp, 5)
five_TRs <- tr_unfav_sp[which(tr_unfav_sp$ID == protein_id_by_number_of_TR(tr_unfav_sp, 5)[[1]]),]
unique(five_TRs$protein_name)
unique(five_TRs$prot_function)

## ------------------------------------------------------------------------
summary(tr_fav_sp$l_effective)
unique(tr_fav_sp$l_type)

## ------------------------------------------------------------------------
summary(tr_unfav_sp$l_effective)
unique(tr_unfav_sp$l_type)

## ------------------------------------------------------------------------
summary(tr_wnt_sp$l_effective)
unique(tr_wnt_sp$l_type)

## ------------------------------------------------------------------------
summary(tr_fav_sp$n_effective)
summary(tr_unfav_sp$n_effective)
summary(tr_wnt_sp$n_effective)

## ------------------------------------------------------------------------
summary(tr_fav_sp$total_repeat_length)
summary(tr_unfav_sp$total_repeat_length)
summary(tr_wnt_sp$total_repeat_length)

## ------------------------------------------------------------------------
tr_wnt_sp$protein_name[which(tr_wnt_sp$total_repeat_length > 15)]
unique(tr_wnt_sp$prot_function[which(tr_wnt_sp$total_repeat_length > 15)])

## ------------------------------------------------------------------------
tr_fav_sp$protein_name[which(tr_fav_sp$total_repeat_length > 15)]
unique(tr_fav_sp$prot_function[which(tr_fav_sp$total_repeat_length > 15)])

## ------------------------------------------------------------------------
tr_unfav_sp$protein_name[which(tr_unfav_sp$total_repeat_length > 15)]
unique(tr_unfav_sp$prot_function[which(tr_unfav_sp$total_repeat_length > 15)])

## ------------------------------------------------------------------------
sel_var <- c("ID", "begin", "msa_original", "repeat_region_length", "l_type", "prot_name", "protein_name", "gene_names", "prot_function")
tr_wnt_sp[which(tr_wnt_sp$ID == "Q9UJU2"),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[which(tr_wnt_sp$ID == "Q9HCS4"),sel_var]
tr_wnt_sp[which(tr_wnt_sp$ID == "P36402"),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("APC", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("GSK", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("KC", tr_wnt_sp$prot_name),sel_var]
tr_fav_sp[grepl("KC", tr_fav_sp$prot_name),sel_var]
tr_unfav_sp[grepl("KC", tr_unfav_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("AXIN", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("DVL", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("FOX", tr_wnt_sp$prot_name),sel_var]
tr_fav_sp[grepl("FOX", tr_fav_sp$prot_name),sel_var]
tr_unfav_sp[grepl("FOX", tr_unfav_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[which(tr_wnt_sp$ID == "O00358"),sel_var]
tr_unfav_sp[which(tr_unfav_sp$ID == "P10070"),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("FZ", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("WNT", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("LRP", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
tr_wnt_sp[grepl("PTK", tr_wnt_sp$prot_name),sel_var]
tr_wnt_sp[grepl("ROR", tr_wnt_sp$prot_name),sel_var]

## ------------------------------------------------------------------------
# url to UniProt/Swiss-Prot querying enzyme commissions for protein kinases
url_kin <- "https://www.uniprot.org/uniprot/?query=ec:2.7.10.-%20OR%20ec:2.7.11.-%20OR%20ec:2.7.12.-%20OR%20ec:2.7.13.-%20OR%20ec:2.7.14.-%20OR%20ec:2.7.99.-&format=fasta&sort=score&fil=proteome:UP000005640%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22"
sp_kinIDs <- load_kinome(url = url_kin, path = dest_file_kin, OnlyIDs = TRUE)

# No. of protein kinases in Human Proteome
length(sp_kinIDs)
length(sp_kinIDs) / nrow(sp_all_fav)
# No. protein kinases in Wnt-Pathway
sum(sp_kinIDs %in% tr_wnt$ID)
sum(sp_kinIDs %in% tr_wnt$ID) / nrow(tr_wnt)
# No. protein kinases in CRC favorable proteins
sum(sp_kinIDs %in% tr_fav$ID)
# No. protein kinases in CRC unfavorable proteins
sum(sp_kinIDs %in% tr_unfav$ID)

(sum(sp_kinIDs %in% tr_fav$ID) + sum(sp_kinIDs %in% tr_unfav$ID))/(nrow(tr_fav)+nrow(tr_unfav))

## ------------------------------------------------------------------------
dest_file_kinext <- "/home/matteo/polybox/MSc_ACLS/master_thesis/TRALResultAnalysis/data/swissprot_human_extkinome.tsv"
url_extkin <- "https://www.uniprot.org/uniprot/?query=ec:2.7.-.-%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
sp_extkinIDS <- load_kinome(url = url_extkin, path = dest_file_kinext, OnlyIDs = TRUE)
# No. of protein kinases in Human Proteome
length(sp_extkinIDS)
length(sp_extkinIDS) / nrow(sp_all_fav)
# No. protein kinases in Wnt-Pathway
sum(sp_extkinIDS %in% tr_wnt$ID)
sum(sp_extkinIDS %in% tr_wnt$ID) / nrow(tr_wnt)
# No. protein kinases in CRC favorable proteins
sum(sp_extkinIDS %in% tr_fav$ID)
# No. protein kinases in CRC unfavorable proteins
sum(sp_extkinIDS %in% tr_unfav$ID)

(sum(sp_extkinIDS %in% tr_fav$ID) + sum(sp_extkinIDS %in% tr_unfav$ID))/(nrow(tr_fav)+nrow(tr_unfav))

## ------------------------------------------------------------------------
top20_unfav_genes <- c("LRCH4", "POFUT2", "CLK3", "EGFL7", "DPP7", "HSH2D", "ASB6", "SPAG4", "EXOC3L4", "HSPA1A", "PAQR6", "FAM69B", "CRACR2B", "ARHGAP4", "NPDC1", "DAPK1", "CNPY3", "ARL8A", "INAFM1", "RHBDD2")
(getTRbyGene(tr_sp = tr_unfav_sp, genes = top20_unfav_genes))

top20_fav_genes <- c("RBM3", "NOL11", "USP53", "TEX2", "HOOK1", "ZYG11B", "HSPA8", "DLAP", "SORT1", "DDX46", "FBXO7", "ABCD3", "NGLY1", "PARS2", "CLCC1", "AP3B1", "PRPSAP1", "PSMA5", "GRSF1", "CD274")
(getTRbyGene(tr_fav_sp, genes = top20_fav_genes))

## ------------------------------------------------------------------------
lina_top6_genes <- c("CDX2", "CLK3", "CNPY3", "CRACR2B", "DPP7", "HSH2D")
lina_top6_prots <- c("Q9NQX5", "H3BVF8", "A0A0C4DFY8", "Q8N4Y2", "R4GMV4", "Q96JZ2")

## ------------------------------------------------------------------------
lina_top6_genes %in% top20_unfav_genes

## ------------------------------------------------------------------------
(getTRbyGene(tr_sp = tr_unfav_sp, genes = lina_top6_genes)[1])

## ------------------------------------------------------------------------
# (aa_freq_TR <- AAfreq_in_TR(Q9UGU0))

## ------------------------------------------------------------------------
# # Not mentioned in Uversky's paper: "B", "O", "U", "Z", "X". These guys might need to fit in with the rest (if possible, as some of them represent multiple aa.)
# aa_order_promoting_to_disorder_promoting = c("C", "W", "I", "Y", "F", "L", "H", "V", "N", "M", "R", "T", "D", "G", "A", "K", "Q", "S", "E", "P", "B", "O", "U", "Z", "X")
# # Sort AA according to their disorder promoting potential
# aa_freq_TR <- aa_freq_TR[match(aa_order_promoting_to_disorder_promoting, aa_freq_TR$aa),]
# colnames(aa_freq_TR) <- c("aa_freq_tr", "aa", "aa_ratio_tr")

## ------------------------------------------------------------------------
# data("AAfreqSP")

## ------------------------------------------------------------------------
# prot_seq <- download_prot_sequence("Q9UGU0")
# aa_freq_prot <- AAfreq_in_prot(prot_seq)
# (aa_freq_prot <- aa_freq_prot[match(aa_order_promoting_to_disorder_promoting, aa_freq_prot$aa),])
# colnames(aa_freq_prot) <- c("aa_freq_sp", "aa", "aa_ratio_sp")

## ------------------------------------------------------------------------
# # we combine the datasets
# aa_freq_prot <- aa_freq_prot[1:20,]
# disorderpropensity <- AAfreqSP$disorderpropensity
# df <- cbind(aa_freq_TR[1:20,], disorderpropensity, aa_freq_prot)
# df <- df[ , !(names(df) %in% c("aa"))]
# 
# p <- ggplot(df, aes(x= aa_ratio_sp, y = aa_ratio_tr, size = disorderpropensity))+
#   geom_point()+
#   labs(x= "AA Background Frequency",
#        y = "AA Frequency in TRs")+
#   guides(size=guide_legend(title="Disorderpropensity"))+
#   theme_minimal()
# p <- beautifier(p, x.axis.text.angle = 0)
# p

