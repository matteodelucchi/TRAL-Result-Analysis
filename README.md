# TRAL Results Analysis

This package provides a reproducible framework of the functional analysis of colorectal cancer associated proteins containing tandem repeats.

* `TRs_favorable_proteins_CRC_sp.tsv`, `TRs_unfavorable_proteins_CRC_sp.tsv`: For each gene associated by the genomeatlas with colorectal cancer the expression data was analysed for tandem repeats with TRAL. 

* `TRs_Wnt_proteins_CRC_sp.tsv`: Each human protein associated with the Wnt pathway in UniProt was analysed with TRAL for tandem repeats. These files contain the results of TRAL.

TRAL results contain i.a. the protein ID, the MSA of the TR region, the begin of the TR region, its length, number of repeats, divergence, p-value, etc.

Further datasets are included which are used for the analysis:

* `swissprot_human.tsv`: Contains information from all human Swiss-Prot proteins. 

* `swissprot_human_kinome.tsv`: Contains information from all human Swiss-Prot protein kinases. 

In the vignette `GeneralOverviewofTandemRepeats.html` all results are reported toghether with an introduction into the topic. Furthermore, its explained how the TRAL results are obtained.

## Installation
To reproduce the results reported in the above mentioned vignette, install the package and adapt the settings in the `.Rmd' version of the vignette according to your settings.

```R
# Install the development version from GitHub
devtools::install_github("matteodelucchi/TRALResultAnalysis")
```

Note, to get the most recent results, download the provided files. To get the exact same results as reported in the vignette, disable the download of any file.
This has to be done manually and will not be automated. To reduce wrong results - [Know the data](http://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/reproducibility-pt-1/)!


