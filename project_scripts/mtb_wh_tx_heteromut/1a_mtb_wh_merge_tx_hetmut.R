#----------------------------------------------------#
# MSc Health Data Science 2023                       #
# Dataset: Mtb Within-host genetic & tx data         #
# Data Linkage: Hetero-resistant strains             #
#----------------------------------------------------#

library(tidyverse)
library(readxl)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
treatment <- read_excel('data/mtb_wh/treatment_table_by_patientID.xlsx')
metadata <- read.csv('data/mtb_wh/mtb_wh_extended.v1.basic_metadata.mi.qc.csv', sep = '\t', header = T)
all_strains <- read.csv('data/mtb_wh/mtb_wh_extended.v1.all_lineages.output_table.mp.csv', sep = '\t', header = T)
gene <- read.table('data/mtb_wh/gene_table_by_strainID.txt', sep = '\t', header = T) # reference table 

# hetero-resistant strains
hetero_mut_raw <- read.csv('data/mtb_wh/mtb_wh_extended.v1.mrca.lofreq_snpEff_annotated_variants.h37rv.in_ref_VCF.filtered.AF0.25_v2.csv', sep = ',', header = T)

# 2. clean variables needed to link datasets----
metadata <- metadata %>% 
  group_by(isolate_id) %>% 
  add_count(isolate_id, name = 'num_isolates') %>% 
  mutate(num_studyid = n_distinct(study_id),
         isolate_n = 1:n()) %>% 
  filter(num_isolates == 1 |  (num_isolates > 1 & num_studyid == 1 & isolate_n == 1) | (num_isolates > 1 & num_studyid == 2 & study_id == 'tb_portals')) %>%
  arrange(strain_id) %>%
  ungroup() %>% 
  select(patient_id, strain_id, geographical_location, study_id) %>% 
  filter(!duplicated(strain_id))

# 3. clean treatment data----
# keep only patient IDs with treatment data
treatment <- treatment %>% 
  filter(regimen_drug != 'NR' & regimen_drug != '{}')

# want to keep only unique patient_id - regimen_drug pairs
treatment <- treatment %>% 
  group_by(isolate_id) %>% 
  mutate(flag = +(n() > 1))

treatment <- treatment[!duplicated(treatment[c("patient_id","regimen_drug")]),]

# then recode treatment vars to include all regimens given to patient 
# code will be different for observations where isolate ID is unique vs where isolate ID is duplicated
unique_isoID <- treatment %>% 
  filter(flag==0)

dup_isoID <- treatment %>% 
  filter(flag==1)

# observations where isolate ID is duplicated
dup_isoID <- dup_isoID %>% 
  select(c(patient_id, regimen_drug, dose, isolate_id)) %>% 
  mutate(regimen_n = 1:n()) %>% 
  pivot_wider(id_cols = c(patient_id, isolate_id),
              names_from = c(regimen_n),
              values_from = c(regimen_drug),
              names_prefix = 'regimen',
              names_sep = '.',
              names_vary = 'slowest')

# re-define treatment variables to include all regimens dispensed over course of condition
drugs <- c('H', 'R', 'S', 'E', 'Z', 'Fq', 'Cm', 'Km', 'Am', 'ART', 'CPT', 'Ofx', 'Lfx', 'Mfx', 'Pto', 'Cs', 'Amx/Clv', 'Mb', 'Dld', 'Bdq', 'Ipm/Cln', 'Lzd', 'Cfz', 'Clr', 'Eto', 'Trd', 'Pas', 'Rifapentine', 'vitamin B6', 'Pcz')

for (drug in drugs) {
  dup_isoID[,drug] <- case_when(grepl(paste0(drug, ','), dup_isoID$regimen1,  fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen2, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen3, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen4, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen5, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen6, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen7, fixed = T) | grepl(paste0(drug, ','), dup_isoID$regimen8, fixed = T) ~ '1',
                                grepl(paste0(drug, '}'), dup_isoID$regimen1,  fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen2, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen3, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen4, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen5, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen6, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen7, fixed = T) | grepl(paste0(drug, '}'), dup_isoID$regimen8, fixed = T) ~ '1',
                                T ~ '0')
}

dup_isoID <- dup_isoID %>% 
  unite(col = 'regimen_drug',
        regimen1:regimen8,
        sep = ';',
        na.rm = T)

# observations where isolate ID is unique
unique_isoID <- unique_isoID %>% 
  select(c(patient_id, regimen_drug, dose, isolate_id)) %>% 
  group_by(patient_id) %>% 
  mutate(regimen_n = 1:n()) %>% 
  pivot_wider(id_cols = c(patient_id),
              names_from = c(regimen_n),
              values_from = c(regimen_drug),
              names_prefix = 'regimen',
              names_sep = '.',
              names_vary = 'slowest') 

# re-define treatment variables to include all regimens dispensed over course of condition
drugs <- c('H', 'R', 'S', 'E', 'Z', 'Fq', 'Cm', 'Km', 'Am', 'ART', 'CPT', 'Ofx', 'Lfx', 'Mfx', 'Pto', 'Cs', 'Amx/Clv', 'Mb', 'Dld', 'Bdq', 'Ipm/Cln', 'Lzd', 'Cfz', 'Clr', 'Eto', 'Trd', 'Pas', 'Rifapentine', 'vitamin B6', 'Pcz')

for (drug in drugs) {
  unique_isoID[,drug] <- case_when(grepl(paste0(drug, ','), unique_isoID$regimen1,  fixed = T) | grepl(paste0(drug, ','), unique_isoID$regimen2, fixed = T) | grepl(paste0(drug, ','), unique_isoID$regimen3, fixed = T) | grepl(paste0(drug, ','), unique_isoID$regimen4, fixed = T) | grepl(paste0(drug, ','), unique_isoID$regimen5, fixed = T) | grepl(paste0(drug, ','), unique_isoID$regimen6, fixed = T) ~ '1',
                                   grepl(paste0(drug, '}'), unique_isoID$regimen1,  fixed = T) | grepl(paste0(drug, '}'), unique_isoID$regimen2, fixed = T) | grepl(paste0(drug, '}'), unique_isoID$regimen3, fixed = T) | grepl(paste0(drug, '}'), unique_isoID$regimen4, fixed = T) | grepl(paste0(drug, '}'), unique_isoID$regimen5, fixed = T) | grepl(paste0(drug, '}'), unique_isoID$regimen6, fixed = T) ~ '1',
                                   T ~ '0')
}

unique_isoID <- unique_isoID %>% 
  unite(col = 'regimen_drug',
        regimen1:regimen6,
        sep = ';',
        na.rm = T)

# combine two sub-datasets into one that is unique at patient_ID
treatment_clean <- rbind(unique_isoID, dup_isoID) %>% 
  select(patient_id, regimen_drug, H:Pcz)

# convert drug vars into numeric vars
treatment_clean[,3:32] <- sapply(treatment_clean[,3:32],as.numeric)

# 4. clean mutation data----
all_strains <- all_strains %>% 
  select(host_id) %>% 
  rename(strain_id = host_id) %>% 
  arrange(strain_id)

hetero_mut <- hetero_mut_raw %>% 
  select(patient_id, isolate_id, locus_tag) %>% 
  group_by(patient_id, isolate_id) %>% 
  mutate(loci_n = 1:n()) %>% 
  pivot_wider(id_cols = c(patient_id, isolate_id),
              names_from = c(loci_n),
              names_prefix = 'locus_tag',
              values_from = c(locus_tag),
              names_sep = '.') %>% 
  unite(col = 'mutated_locus',
        locus_tag1:locus_tag246,
        sep = '',
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
        isolate1:isolate49,
        sep = '',
        na.rm = T) %>% 
  ungroup()

# need to create binary variables for each locus on Mtb genome
all_loci <- names(gene[,3:4189])

for (loci in all_loci) {
  hetero_mut[,loci] <- case_when(grepl(loci, hetero_mut$mutated_locus, fixed = T) ~ 1,
                                    T ~ 0)
}

hetero_mut <- hetero_mut %>% 
  filter(mutated_locus != '') %>% 
  rename(strain_id = patient_id)

# merge list of all strains into gene data
hetero_mut <- hetero_mut %>% 
  merge(all_strains, by = 'strain_id', all.x = T, all.y = T) %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  arrange(strain_id)

# 5. merge treatment data with strain ids----
treatment_clean <- treatment_clean %>% 
  merge(metadata, by = 'patient_id') %>% 
  select(c(patient_id, strain_id, geographical_location, study_id, regimen_drug, H:Pcz))

# 6. merge treatment data with gene/mutation data----
treatment_hetero_mutation <- treatment_clean %>% 
  merge(hetero_mut, by = 'strain_id') %>% 
  select(c(patient_id, strain_id, geographical_location, study_id, regimen_drug:Pcz, mutated_locus:Rv3924c)) %>% 
  arrange(patient_id)

# 7. recode Fq and Eto/Pto----
treatment_hetero_mutation <- treatment_hetero_mutation %>% 
  mutate(Fq_recode = case_when(Fq == 1 ~ 1, # recode Fq to include Ofx and Lfx
                               Ofx == 1 ~ 1,
                               Lfx == 1 ~ 1,
                               Mfx == 1 ~ 1,
                               T ~ 0),
         Eto_recode = case_when(Eto == 1 ~ 1,
                                Pto == 1 ~ 1,
                                T ~ 0)) %>% 
  select(patient_id, strain_id, geographical_location, study_id,
         regimen_drug:Z, Fq_recode, Cm:Am, Cs:Clr, Eto_recode, Trd:Pas, 
         mutated_locus:Rv3924c) %>% 
  rename(Amx = 'Amx/Clv',
         Ipm = 'Ipm/Cln',
         Fq = Fq_recode,
         Eto = Eto_recode)

# 7. save dataset----
saveRDS(treatment_hetero_mutation, file = "~/Desktop/LSHTM/Summer Project/analysis/datasets/tx_hetmut.Rds")
