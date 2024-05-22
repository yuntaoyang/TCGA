# TCGA
This repository contains scripts for downloading TCGA datasets and performing survival analysis.
1. Download TCGA datasets
* Check `TCGA_Download.R`.
* Output data
    * Raw reads count `./data/TCGA-BLCA_counts.csv`.
    * Clinical data `./data/TCGA-BLCA_clinical.csv`.
    * Biospecimen data `./data/TCGA-BLCA_biospecimen.csv`.
2. Perform survival analysis using maxstat
* Check `TCGA_Survival.R`.
* Input data
    * Raw reads count `./data/TCGA-BLCA_counts.csv`.
    * Clinical data `./data/TCGA-BLCA_clinical.csv`.
    * Biospecimen data `./data/TCGA-BLCA_biospecimen.csv`.
* Output data
    * Survival data `./out/TCGA-BLCA_survival.csv`.
    * Kaplan-Meier plot `./out/TCGA-BLCA_survival.pdf`.