## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(123)

## ----sessioninfo, echo=FALSE---------------------------------------------
sessionInfo()

## ----Housekeeping, message=FALSE, warning=FALSE, collapse=TRUE-----------
library(TRALResultAnalysis)
tr_crcfavorable_path <- "./inst/extdata/TRs_favorable_proteins_CRC_sp.tsv"
tr_crcunfavorable_path <- "./inst/extdata/TRs_unfavorable_proteins_CRC_sp.tsv"
tr_wnt_path <- "./inst/extdata/TRs_Wnt_proteins_CRC_sp.tsv"
dest_file <-"./data/swissprot_human.tsv"

