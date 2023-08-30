# LSHTM MSc Health Data Science Thesis
This repository contains the project source code for the LSTHM MSc Health Data Science thesis titled "*Investigating anti-tuberculosis drug regimens and Mycobacterium tuberculosis strain characteristics as drivers of drug resistance acquisition during treatment*."

### data sources
This project utilized data from TB Portals (https://tbportals.niaid.nih.gov/) and publicly available datasets from published studies.

### source code
The R scripts can be found in the folder titled "project_scripts" and are divided into two sub-folders, each pertaining to either the primary or secondary analysis of the project. A description of each folder and its contents can be found below.

#### primary analysis: mtb_wh_fixedmut
This folder contains the R scripts for the generation of the primary analysis dataset and the primary analysis of the project, where the outcome was the acquisition of fixed mutations in the *Mtb* genome that confer drug resistance. Each R script is labeled in the order that it should be run: scripts labeled with a prefix starting in 1 are used to clean and merge input datasets, generate new variables, and generate final analysis datasets; scripts labeled with a prefix starting in 2 are used to run project analyses.

* **1a_mtb_wh_gen_mutationtable**: generates the input file that determines which *Mtb* strains had a fixed *de novo* mutation (i.e., originated during treatment) in each gene in the *Mtb* genome
* **1b_mtb_wh_merge_tx_fixedmut**: generates a dataset that contains all the *Mtb* strains for which both treatment and fixed mutation data were available.
* **1c_mtb_wh_merge_tx_fixedmut_tbprof**: uses the dataset generated in 1b and merges *TBProfiler* data for all strains and identifies each strain's baseline drug resistance profile.
* **1d_mtb_wh_gen_tx_fixedmut_neffectdrugs**: uses the dataset generated in 1d and generates a variable measuring the total number of drugs used during treatment and the number of effective drugs used during treatment.

* **2a_mtb_wh_datasummary**: generates summary statistics for the main variables used in the analysis; also generates figures included in the final thesis report.
* **2b_mtb_wh_freqtable_tx_fixedmut**: generates the output files that were used to create a frequency table measuring the frequency of drug resistance acquisition for individual anti-tuberculosis drugs.
* **2c_mtb_wh_reg_tx_fixedmut**: conducts the analysis measuring the association between the use of each anti-tuberculosis drug during treatment and likelihood of acquiring fixed mutation on *gyrA*, a gene known to confer resistance to fluoroquinolones.
* **2d_mtb_wh_reg_neffectdrug_fixedmut_v2**: conducts the analysis measuring the associations between the number of effective drugs used during treatment, strain lineage, baseline drug resistance, and total number of drugs used during treatment and likelihood of acquiring drug resistance.
* **2e_mtb_wh_reg_fixedmut_txtime**: conducts the analysis measuring the assocations between the number of effective drugs used during treatment and total treatment duration and the likelihood of acquiring drug resistance.

#### secondary analysis: mtb_wh_heteromut
This folder contains the R scripts for the generation of the secondary analysis dataset and the secondary analysis of the project, where the outcome was the acquisition of heterozygous mutations in the *Mtb* genome that confer drug resistance. Each R script is labeled in the order that it should be run: scripts labeled with a prefix starting in 1 are used to clean and merge input datasets, generate new variables, and generate final analysis datasets; scripts labeled with a prefix starting in 2 are used to run project analyses.

* **1a_mtb_wh_merge_tx_hetmut**: generates a dataset that contains all the *Mtb* strains for which both treatment and heterozygous mutation data were available.
* **1b_mtb_wh_merge_tx_hetmut_tbprof**: uses the dataset generated in 1b and merges *TBProfiler* data for all strains and identifies each strain's baseline drug resistance profile.
* **1c_mtb_wh_gen_tx_hetmut_neffectdrugs**: uses the dataset generated in 1d and generates a variable measuring the total number of drugs used during treatment and the number of effective drugs used during treatment.

* **2a_mtb_wh_freqtable_tx_hetmut**: generates the output files that were used to create a frequency table measuring the frequency of drug resistance acquisition for individual anti-tuberculosis drugs.
* **2b_mtb_wh_reg_tx_hetmut**: conducts the analysis measuring the association between the use of each anti-tuberculosis drug during treatment and likelihood of acquiring heterozygous mutation on *gyrA*, a gene known to confer resistance to fluoroquinolones.
* **2c_mtb_wh_reg_neffectdrug_hetmut_v2**: conducts the analysis measuring the associations between the number of effective drugs used during treatment, strain lineage, baseline drug resistance, and total number of drugs used during treatment and likelihood of acquiring drug resistance.
* **2d_mtb_wh_reg_hetmut_txtime**: conducts the analysis measuring the assocations between the number of effective drugs used during treatment and total treatment duration and the likelihood of acquiring drug resistance.
