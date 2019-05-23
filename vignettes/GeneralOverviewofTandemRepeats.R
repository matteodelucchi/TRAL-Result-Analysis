## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(123)

## ----Housekeeping, message=FALSE, warning=FALSE, collapse=TRUE-----------
library(TRALResultAnalysis)
tr_crcfavorable_path <- "/home/matteo/polybox/MSc_ACLS/master_thesis/TRALResultAnalysis/inst/extdata/TRs_favorable_proteins_CRC.tsv"
tr_crcunfavorable_path <- "/home/matteo/polybox/MSc_ACLS/master_thesis/TRALResultAnalysis/inst/extdata/TRs_unfavorable_proteins_CRC.tsv"
dest_file <-"/home/matteo/polybox/MSc_ACLS/master_thesis/TRALResultAnalysis/data/swissprot_human.tsv"

## ----eval=TRUE, echo=FALSE, message=TRUE, warning=TRUE-------------------
tr_fav = load_tr_annotations(tr_crcfavorable_path)
tr_unfav = load_tr_annotations(tr_crcunfavorable_path)

## ----overview favorable, collapse=TRUE-----------------------------------
str(tr_fav)
head(tr_fav)

## ----overview unfavorable------------------------------------------------
str(tr_unfav)

## ----Download sp_all if it doesnt exist----------------------------------
sp_url <- "https://www.uniprot.org/uniprot/?query=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20reviewed:yes&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,virus%20hosts,encodedon,database(Pfam),interactor,comment(ABSORPTION),feature(ACTIVE%20SITE),comment(ACTIVITY%20REGULATION),feature(BINDING%20SITE),feature(CALCIUM%20BIND),comment(CATALYTIC%20ACTIVITY),comment(COFACTOR),feature(DNA%20BINDING),ec,comment(FUNCTION),comment(KINETICS),feature(METAL%20BINDING),feature(NP%20BIND),comment(PATHWAY),comment(PH%20DEPENDENCE),comment(REDOX%20POTENTIAL),rhea-id,feature(SITE),comment(TEMPERATURE%20DEPENDENCE)&sort=organism"

if(!file.exists(dest_file)){
    download.file(sp_url, destfile = dest_file)
}

## ----load sp_all---------------------------------------------------------
sp_all_fav <- load_swissprot(dest_file, tr_fav)
sp_all_unfav <- load_swissprot(dest_file, tr_unfav)

## ----add meta data-------------------------------------------------------
tr_unfav_sp = merge(x = tr_unfav, y = sp_all_unfav, by = "ID", all.x = TRUE)
tr_fav_sp = merge(x = tr_fav, y = sp_all_fav, by = "ID", all.x = TRUE)

## ------------------------------------------------------------------------
str(tr_unfav_sp)

## ------------------------------------------------------------------------
length(tr_fav$ID)

## ------------------------------------------------------------------------
length(unique(tr_fav$ID))

## ------------------------------------------------------------------------
table(table(tr_fav$ID))

## ------------------------------------------------------------------------
protein_id_by_number_of_TR(tr_fav, 6)
protein_id_by_number_of_TR(tr_fav, 3)

## ------------------------------------------------------------------------
Q9UGU0 <- tr_fav_sp[which(tr_fav_sp$ID == "Q9UGU0"),]

## ------------------------------------------------------------------------
Q9UGU0$prot_function[1]

## ------------------------------------------------------------------------
Q9UGU0$l_type

## ------------------------------------------------------------------------
TR_location(tr_all_sp = Q9UGU0, byTRtype = FALSE)

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
(aa_freq_TR <- AAfreq_in_TR(Q9UGU0))

## ------------------------------------------------------------------------
# Not mentioned in Uversky's paper: "B", "O", "U", "Z", "X". These guys might need to fit in with the rest (if possible, as some of them represent multiple aa.)
aa_order_promoting_to_disorder_promoting = c("C", "W", "I", "Y", "F", "L", "H", "V", "N", "M", "R", "T", "D", "G", "A", "K", "Q", "S", "E", "P", "B", "O", "U", "Z", "X")
# Sort AA according to their disorder promoting potential
aa_freq_TR <- aa_freq_TR[match(aa_order_promoting_to_disorder_promoting, aa_freq_TR$aa),]
colnames(aa_freq_TR) <- c("aa_freq_tr", "aa", "aa_ratio_tr")

## ------------------------------------------------------------------------
data("AAfreqSP")

## ------------------------------------------------------------------------
prot_seq <- download_prot_sequence("Q9UGU0")
aa_freq_prot <- AAfreq_in_prot(prot_seq)
(aa_freq_prot <- aa_freq_prot[match(aa_order_promoting_to_disorder_promoting, aa_freq_prot$aa),])
colnames(aa_freq_prot) <- c("aa_freq_sp", "aa", "aa_ratio_sp")

## ------------------------------------------------------------------------
# we combine the datasets
aa_freq_prot <- aa_freq_prot[1:20,]
disorderpropensity <- AAfreqSP$disorderpropensity
df <- cbind(aa_freq_TR[1:20,], disorderpropensity, aa_freq_prot)
df <- df[ , !(names(df) %in% c("aa"))]

p <- ggplot(df, aes(x= aa_ratio_sp, y = aa_ratio_tr, size = disorderpropensity))+
  geom_point()+
  labs(x= "AA Background Frequency",
       y = "AA Frequency in TRs")+
  guides(size=guide_legend(title="Disorderpropensity"))+
  theme_minimal()
p <- beautifier(p, x.axis.text.angle = 0)
p

