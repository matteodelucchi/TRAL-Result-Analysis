library(RCurl)
library(Biostrings)
library(tidyr)

#################
### Download all protein sequences from swissprot and store them as fasta
#################
#sp_url <- "https://www.uniprot.org/uniprot/?query=reviewed:yes&format=fasta&limit=10&columns=id,entry%20name,protein%20names,sequence" #limited to the first 10 entries. Use for debugging.
sp_url <- "https://www.uniprot.org/uniprot/?query=reviewed:yes&format=fasta&columns=id,entry%20name,protein%20names,sequence"
aa_sp <- getURL(sp_url)
temp <- tempfile() # store it in a temporary file
write.table(aa_sp, file = temp, sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)
aa_SP <- readAAStringSet(temp, format = "fasta")
unlink(temp)

#################
### Count all AAs in all swissprots
#################
AAfreqSP <- AAfreq_in_prot(aa_SP)
colnames(AAfreqSP) <- c("aa_freq_sp", "aa", "aa_ratio_sp")

#################
### Add disorder information
#################
# Uversky's paper: http://www.tandfonline.com/doi/full/10.4161/idp.24684
# Not mentioned in Uversky's paper: "B", "O", "U", "Z", "X". These guys might need to fit in with the rest (if possible, as some of them represent multiple aa.)
aa_order_promoting_to_disorder_promoting = c("C", "W", "I", "Y", "F", "L", "H", "V", "N", "M", "R", "T", "D", "G", "A", "K", "Q", "S", "E", "P", "B", "O", "U", "Z", "X")
# Sort AA according to their disorder promoting potential
AAfreqSP <- AAfreqSP[match(aa_order_promoting_to_disorder_promoting, AAfreqSP$aa),]
# Add the disorder propensity from Uversky paper
AAfreqSP$disorderpropensity <- c(0.00, 0.004, 0.090, 0.113, 0.117, 0.195, 0.259, 0.263, 0.285, 0.291, 0.394, 0.401, 0.407, 0.437, 0.450, 0.588,0.665,0.713, 0.781,1.000, NA, NA, NA, NA, NA)
# calculate disorder propensity manually
#TODO
# remove all AA for which we don't have disorderpropensity values
AAfreqSP <- AAfreqSP %>%
  drop_na()

#################
### Save to RData file
#################
# save it as .tsv
write.table(AAfreqSP, file = "./data-raw/AAfreqSP.tsv", sep = "\t", row.names = FALSE)
# make the data beeing delivered with the package
usethis::use_data(AAfreqSP, overwrite=TRUE)
