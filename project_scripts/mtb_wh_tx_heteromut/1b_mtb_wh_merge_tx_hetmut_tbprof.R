#-------------------------------------------------------------#
# MSc Health Data Science 2023                                #
# Dataset: Mtb Within-host genetic & tx data                  #
# Merge in baseline resistance data                           #
#-------------------------------------------------------------#

library(tidyverse)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_hetmut.Rds') 
metadata <- read.csv('data/mtb_wh/mtb_wh_extended.v1.basic_metadata.mi.qc.csv', sep = '\t', header = T)
tbprofiler <- read.table('data/mtb_wh/TBProfiler_resistance_table.txt', sep = '\t', header = T)

tbportals_gen <- read.csv('data/TBPortals/genomic_data/TB_Portals_Published_Genomics_Data_January_2023/TB_Portals_Genomics_January_2023.csv')
tbportals_pat <- read.csv('data/TBPortals/clinical_data/TB_Portals_Published_Clinical_Data_January_2023/TB_Portals_Patient_Cases_January_2023.csv')

# use tbportals to get patient IDs and isolate IDs for those in tbportals
tbportals <- tbportals_pat %>% 
  merge(y = tbportals_gen, by = 'condition_id') %>% 
  select(patient_id, condition_id, sra_id) %>% 
  arrange(patient_id, sra_id) %>% 
  rename(isolate_id = sra_id)

# then filter tbportals data for those in master dataset
tbportals <- tbportals %>% 
  merge(y = data, by = 'patient_id') %>% 
  select(patient_id, strain_id, isolate_id) %>% 
  rename(sample = 'isolate_id')

# then get strain IDs for non tbportals data
incl_strain_id <- data %>% 
  filter(study_id != 'tb_portals') %>% 
  merge(y = metadata, by = c('patient_id', 'strain_id', 'geographical_location', 'study_id')) %>% 
  select(patient_id, strain_id, isolate_id) %>%
  filter(!duplicated(isolate_id)) %>% # some isolates are listed 2x in metadata but have same data
  rename(sample = isolate_id)

# 2. clean geography data----
# recode geography values
table(data$geographical_location)

data <- data %>%
  rename(geography = geographical_location) %>% 
  mutate(geography_recode = case_when(geography == 'Abkhazia, Republic of Georgia' ~ 'Georgia',
                                      geography == 'Beijing, China' | geography == 'China: shenzhen' | geography == 'China:Jiangsu' | geography == 'China:Sichuan' | geography == 'China:Zhejiang' | geography == 'Zhejiang Province, China' ~ 'China',
                                      geography == 'British Columbia, Canada' | geography == 'Canada: Ontario' ~ 'Canada',
                                      geography == 'Cape_Town, South Africa' | geography == 'Francistown, South Africa' | geography == 'Harare, South Africa' | geography == 'Johannesburg, South Africa' | geography == 'Khayelitsha, Cape Town, South Africa' | geography == 'Marondera, South Africa' | geography == 'South Africa: Durban' | geography == 'South Africa: Hlabisa' | geography == 'South Africa: Kwazulu-Natal' | geography == 'Western Cape, South Africa' ~ 'South Africa',
                                      geography == 'India: Chennai' ~ 'India',
                                      geography == 'Japan:Okinawa' | geography == 'Japan:Osaka' | geography == 'Japan:Osaka, Osaka' | geography == 'Japan:Tokyo' | geography == 'Japan:Yamagata' | geography == 'Japan:Hyogo, Kobe' ~ 'Japan',
                                      geography == 'Kampala, Uganda' ~ 'Uganda',
                                      geography == 'Karakalpakstan, Uzbekistan' ~ 'Uzbekistan',
                                      geography == 'Karonga, Malawi' ~ 'Malawi',
                                      geography == 'London, UK' | geography == 'Midlands, UK' ~ 'UK',
                                      geography == 'Minsk, Belarus' ~ 'Belarus',
                                      geography == 'Russia: Kaliningrad' ~ 'Russia',
                                      geography == 'Sweden: Stockholm' ~ 'Sweden',
                                      geography == 'Valencia, Spain' ~ 'Spain',
                                      T ~ geography
                                      )) %>% 
  select(patient_id:geography, geography_recode, study_id:Rv3924c)

# 3. clean strain and baseline resistance profile data (TBProfiler)----
# merge in patient_id/strain_id from tbportals/metadata and keep only patients/strains in master dataset
# tbportals
tbportals_tbprof <- tbportals %>% 
  merge(y = tbprofiler, by = 'sample') %>% 
  select(patient_id, strain_id, sample, tb_profiler_strain:resistance_type)
  
# non tbportals
notbportals_tbprof <- incl_strain_id %>% 
  merge(y = tbprofiler, by = 'sample') %>% 
  select(patient_id, strain_id, sample, tb_profiler_strain:resistance_type)

# append two datasets together
tbprofiler_merge <- tbportals_tbprof %>% 
  rbind(notbportals_tbprof) %>% 
  arrange(patient_id, strain_id, sample) %>% 
  # then proceed with process for filtering for baseline resistance profile
  add_count(strain_id, name = 'num_isolates') %>%
  group_by(strain_id) %>%
  mutate(isolate_n = 1:n(),
         resistance_type_n = n_distinct(resistance_type))

# create new data frame with baseline resistance profile 
tbprofiler_baseline <- tbprofiler_merge %>% 
  filter(resistance_type_n == 1) %>% # only include isolates that don't change resistance type over time first
  separate(col = resistant_drug, 
           into = c('drug1', 'drug2', 'drug3', 'drug4', 'drug5', 'drug6', 'drug7', 'drug8', 'drug9', 'drug10'),
           sep = ';',
           remove = F)

# help count number of drugs each isolate is resistant to
count_na_func <- function(x) sum(is.na(x))

tbprofiler_baseline$count_na <- apply(tbprofiler_baseline[,7:16], 1, count_na_func)
tbprofiler_baseline$num_drug_resist <- 10 - tbprofiler_baseline$count_na

tbprofiler_baseline <- tbprofiler_baseline %>% 
  filter((n_distinct(num_drug_resist) == 1 & isolate_n == min(isolate_n)) | 
           (n_distinct(num_drug_resist) > 1 & num_drug_resist == min(num_drug_resist))) %>% 
  add_count(strain_id, name = 'n_samples') %>% 
  filter(n_samples == 1 | n_samples > 1 & isolate_n == min(isolate_n)) %>% 
  select(strain_id:resistant_drug, mutations:resistance_type)

# separate isolates that have different resistance profiles over time
tbprofiler_multi_resist <- tbprofiler_merge %>% 
  filter(resistance_type_n > 1) %>% # include isolates with different resistance profiles
  separate(col = resistant_drug, 
           into = c('drug1', 'drug2', 'drug3', 'drug4', 'drug5', 'drug6', 'drug7', 'drug8', 'drug9', 'drug10'),
           sep = ";",
           remove = F)

tbprofiler_multi_resist$count_na <- apply(tbprofiler_multi_resist[,7:16], 1, count_na_func)
tbprofiler_multi_resist$num_drug_resist <- 10 - tbprofiler_multi_resist$count_na

tbprofiler_multi_resist <- tbprofiler_multi_resist %>% 
  group_by(strain_id) %>% 
  mutate(baseline_profile = min(num_drug_resist)) %>% # determine least number of drugs each isolate is resistant to
  filter(num_drug_resist == baseline_profile) %>% # filter those where the number of drugs is equal to minimum
  mutate(flag1 = case_when(n_distinct(sample) > 1 ~ 1, # flag those where there are still multiple isolates per strain
                           T ~ 0),
         flag2 = case_when(n_distinct(resistance_type) > 1 ~ 1, # flag those where resistance type still changes over time and resistant to same number of drugs
                           T ~ 0),
         keep = case_when(#flag1==0 ~ 1,
                          flag1==0 & flag2==0 & isolate_n==min(isolate_n) ~ 1,
                          flag1==1 & flag2==0 & isolate_n==min(isolate_n) ~ 1,
                          flag1==1 & flag2==1 & (sample=='SRR6807681' | sample=='SRR13232497' | sample=='SRR6807756') ~ 1,
                          T ~ 0)) %>% # flag the isolate kept for each strain
  filter(keep==1) %>% # keep baseline profile for each strain
  select(strain_id:resistant_drug, mutations:resistance_type)

# append to baseline data frame
tbprofiler_baseline <- tbprofiler_baseline %>% 
  rbind(tbprofiler_multi_resist) %>% 
  arrange(strain_id) %>% 
  ungroup()

# gen main_lineage var (mixed strains will have value 'mixed')
tbprofiler_baseline <- tbprofiler_baseline %>% 
  mutate(main_lineage = case_when(mixed_strain == 'yes' ~ 'mixed',
                                  tb_profiler_strain == 'lineage1' | grepl('lineage1.', tb_profiler_strain, fixed = T)  ~ 'lineage1',
                                  tb_profiler_strain == 'lineage2' | grepl('lineage2.', tb_profiler_strain, fixed = T)  ~ 'lineage2',
                                  tb_profiler_strain == 'lineage3' | grepl('lineage3.', tb_profiler_strain, fixed = T)  ~ 'lineage3',
                                  tb_profiler_strain == 'lineage4' | grepl('lineage4.', tb_profiler_strain,  fixed = T) ~ 'lineage4',
                                  tb_profiler_strain == 'lineage5' | grepl('lineage5.', tb_profiler_strain,  fixed = T) ~ 'lineage5',
                                  tb_profiler_strain == 'lineage6' | grepl('lineage6.', tb_profiler_strain,  fixed = T) ~ 'lineage6',
                                  tb_profiler_strain == 'lineage9' | grepl('lineage9.', tb_profiler_strain,  fixed = T) ~ 'lineage9',
  ))

# 4. merge tbprofiler data with main dataset----
data <- data %>% 
  merge(y = tbprofiler_baseline, by = 'strain_id') %>% 
  select(strain_id, patient_id, study_id, sample, geography, geography_recode, 
         tb_profiler_strain, main_lineage, mixed_strain, resistant_drug, mutations, resistance_type,
         regimen_drug:Rv3924c)

# 5. save dataset----
saveRDS(data, file = "~/Desktop/LSHTM/Summer Project/analysis/datasets/tx_hetmut_tbprof.Rds")
