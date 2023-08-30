#----------------------------------------------------#
# MSc Health Data Science 2023                       #
# Dataset: Mtb Within-host genetic & tx data         #
# Data validation: recreate gene_table_by_strainID   #
#----------------------------------------------------#

library(tidyverse)
library(readxl)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
metadata <- read.csv('data/mtb_wh/mtb_wh_extended.v1.basic_metadata.mi.qc.csv', sep = '\t', header = T)
gene <- read.table('data/mtb_wh/gene_table_by_strainID.txt', sep = '\t', header = T) # reference table 

# mutation data (fixed mutations)
fixed_mut_raw <- read.csv('data/mtb_wh/filtered_mutations_added_columns.csv', sep = ',', header = T)

# 1. clean mutation data----
fixed_mut <- fixed_mut_raw %>%  
  filter(locus_tag != '' & locus_tag != '-') %>% 
  select(patient_id, isolate_id, locus_tag) %>% 
  group_by(patient_id, isolate_id) %>% 
  mutate(loci_n = 1:n()) %>% 
  pivot_wider(id_cols = c(patient_id, isolate_id),
              names_from = c(loci_n),
              names_prefix = 'locus_tag',
              values_from = c(locus_tag),
              names_sep = '.') %>% 
  unite(col = 'mutated_locus',
        locus_tag1:locus_tag130,
        sep = ',',
        na.rm = T) %>%
  arrange(patient_id) %>%
  group_by(patient_id) %>%
  mutate(isolate_n = 1:n()) %>%
  pivot_wider(id_cols = patient_id,
              names_from = isolate_n,
              names_prefix = 'isolate',
              values_from = mutated_locus,
              names_sep = '.') %>%
  unite(col = 'mutated_locus',
        isolate1:isolate26,
        sep = ',',
        na.rm = T) %>%
  ungroup()

# need to create binary variables for each locus on Mtb genome
all_loci <- names(gene[,3:4189])

for (loci in all_loci) {
  fixed_mut[,loci] <- case_when(grepl(loci, fixed_mut$mutated_locus, fixed = T) ~ 1,
                                   T ~ 0)
}

# 2. save dataset----
fixed_mut <- fixed_mut %>% 
  rename(strain_id = patient_id)

saveRDS(fixed_mut, file = "~/Desktop/LSHTM/Summer Project/analysis/datasets/fixed_mutations_nmj.Rds")
