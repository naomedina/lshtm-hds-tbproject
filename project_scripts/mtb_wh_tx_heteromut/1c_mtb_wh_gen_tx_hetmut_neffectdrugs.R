#-------------------------------------------------------------#
# MSc Health Data Science 2023                                #
# Dataset: Mtb Within-host genetic & tx data                  #
# Determine num of effective drugs for each strain            #
#-------------------------------------------------------------#

library(tidyverse)
library(Dict)
library(patchwork)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_hetmut_tbprof.Rds')

# 2. factor resistance type----
data <- data %>% 
  rename(baseline_resistance_type = resistance_type)

data$baseline_resistance_type <- factor(data$baseline_resistance_type,
                                        levels = c('pan-susceptible', 'mono-resistant', 'resistant', 'MDR', 'XDR'))

# tab baseline resistance categories
data %>% 
  group_by(baseline_resistance_type) %>% 
  summarise(n = n())

# pan-susceptible: no predicted drug resistance
# mono-resistant: predicted resistance to one drug
# resistant: predicted resistance to more than one drug, but cannot be classified as MDR
# MDR: predicted resistance to (at least) isoniazid and rifampicin
# XDR: predicted resistance to isoniazid, rifampicin, a fluoroquinolone, and at least one second-line injectable

# 3. create count for number of resistant drugs and determine effective drugs----
# 22 drugs
drugs <- c('Isoniazid', 'Rifampicin', 'Streptomycin', 'Ethambutol', 'Pyrazinamide', 'Fluoroquinolones', 'Capreomycin', 'Kanamycin', 'Amikacin', 'Cycloserine', 'Amoxicillin-clavulanate', 'Delamanid', 'Bedaquiline', 'Imipenem-cilastatin', 'Linezolid', 'Clofazimine', 'Clarithromycin', 'Ethionamide', 'Terizidone', 'Para-aminosalicylic_acid')

data <- data %>% 
  mutate(num_resistant_drugs = 0,
         effective_drugs = '')

for (row in 1:nrow(data)) {
  num_resist = 0
  effective_drugs = c()
  
  for (drug in drugs) {
    
    if (drug == 'Ethionamide') {
      num_resist = case_when(grepl(drug, data[row, 'resistant_drug'], fixed = T) | grepl('Prothionamide', data[row, 'resistant_drug'], fixed = T) ~ num_resist + 1, 
                             T  ~ num_resist)
      if (!grepl(drug, data[row, 'resistant_drug'], fixed = T) & !grepl('Prothionamide', data[row, 'resistant_drug'], fixed = T)) {
        effective_drugs <- append(effective_drugs, drug)
      }
    }
    
    else {
    num_resist = case_when(grepl(drug, data[row, 'resistant_drug'], fixed = T) ~ num_resist + 1, 
                           T  ~ num_resist)
      if (!grepl(drug, data[row, 'resistant_drug'], fixed = T)) {
        effective_drugs <- append(effective_drugs, drug)
      }
    }
    
  data[row, 'num_resistant_drugs'] = num_resist
  data[row, 'effective_drugs'] = paste0(effective_drugs, collapse = ', ')
  
  }
}

# 4. create num of effective drugs var for each strain----
data <- data %>% 
  mutate(num_effective_drugs = length(drugs) - num_resistant_drugs)

# 5. determine number of effective drugs each strain is treated with----
# create a dictionary
tb_drugs <- Dict$new(
  'H' = 'Isoniazid',
  'R' = 'Rifampicin',
  'S' = 'Streptomycin',
  'E' = 'Ethambutol',
  'Z' = 'Pyrazinamide',
  'Fq' = 'Fluoroquinolones',
  'Cm' = 'Capreomycin',
  'Km' = 'Kanamycin',
  'Am' = 'Amikacin',
  'Cs' = 'Cycloserine',
  'Amx' = 'Amoxicillin-clavulanate',
  'Mb' = 'Mycobutin',
  'Dld' = 'Delamanid',
  'Bdq' = 'Bedaquiline',
  'Ipm' = 'Imipenem-cilastatin',
  'Lzd' = 'Linezolid',
  'Cfz' = 'Clofazimine',
  'Clr' = 'Clarithromycin',
  'Eto' = 'Ethionamide',
  'Trd' = 'Terizidone',
  'Pas' = 'Para-aminosalicylic_acid'
)

data <- data %>% 
  mutate(tx_effective_drugs = 0)

for (row in 1:nrow(data)) {
  tx_effective_drugs = 0
  for (var in 14:34) {
    varname <- colnames(data)[var]
    val <- tb_drugs$get(varname)
    tx_effective_drugs = case_when(grepl(val, data[row, 'effective_drugs'], fixed = T) & data[row, varname] == 1 ~ tx_effective_drugs + 1, 
                                   T  ~ tx_effective_drugs)
  }
  data[row, 'tx_effective_drugs'] = tx_effective_drugs
}

# # summarise new var
# data %>% 
#   summarise(min = min(tx_effective_drugs),
#             p25 = quantile(tx_effective_drugs, 0.25),
#             median = median(tx_effective_drugs),
#             mean = mean(tx_effective_drugs),
#             p75 = quantile(tx_effective_drugs, 0.75),
#             max = max(tx_effective_drugs))
# 
# data %>% 
#   group_by(baseline_resistance_type) %>% 
#   summarise(min = min(tx_effective_drugs),
#             p25 = quantile(tx_effective_drugs, 0.25),
#             median = median(tx_effective_drugs),
#             mean = mean(tx_effective_drugs),
#             p75 = quantile(tx_effective_drugs, 0.75),
#             max = max(tx_effective_drugs))
# 
# hist(data$tx_effective_drugs, 
#      col = 'cadetblue3',
#      main = 'Number of effective drugs in treatment regimen',
#      xlab = '')
# 
# data %>% 
#   ggplot(aes(x = baseline_resistance_type, y = tx_effective_drugs)) +
#   geom_boxplot(fill = 'cadetblue3') +
#   xlab('baseline resistance') +
#   ylab('number of effective drugs in treatment regimen') + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "lightgrey"),
#         panel.border = element_blank())

data %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv0006==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv2043c==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv0667==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv3795==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv3423c==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv1908c==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv2780==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv3854c==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(Rv0678==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

data %>% 
  filter(MTB000019==1) %>% 
  group_by(tx_effective_drugs) %>% 
  summarise(n = n())

# 6. determine number of drugs resistant at baseline and used during tx----
data <- data %>% 
  mutate(tx_resistant_drugs = 0)

for (row in 1:nrow(data)) {
  tx_resistant_drugs = 0
  for (var in 14:34) {
    varname <- colnames(data)[var]
    val <- tb_drugs$get(varname)
    tx_resistant_drugs = case_when(grepl(val, data[row, 'resistant_drug'], fixed = T) & data[row, varname] == 1 ~ tx_resistant_drugs + 1, 
                                   T  ~ tx_resistant_drugs)
  }
  data[row, 'tx_resistant_drugs'] = tx_resistant_drugs
}

# # summarise new var
# data %>% 
#   summarise(min = min(tx_resistant_drugs),
#             p25 = quantile(tx_resistant_drugs, 0.25),
#             median = median(tx_resistant_drugs),
#             mean = mean(tx_resistant_drugs),
#             p75 = quantile(tx_resistant_drugs, 0.75),
#             max = max(tx_resistant_drugs))
# 
# data %>% 
#   group_by(baseline_resistance_type) %>% 
#   summarise(min = min(tx_resistant_drugs),
#             p25 = quantile(tx_resistant_drugs, 0.25),
#             median = median(tx_resistant_drugs),
#             mean = mean(tx_resistant_drugs),
#             p75 = quantile(tx_resistant_drugs, 0.75),
#             max = max(tx_resistant_drugs))
# 
# hist(data$tx_resistant_drugs, 
#      col = 'cadetblue3',
#      main = 'Number of resistant drugs in treatment regimen',
#      xlab = '')
# 
# data %>% ggplot(aes(x = baseline_resistance_type, y = tx_resistant_drugs)) +
#   geom_boxplot(fill = 'cadetblue3') +
#   xlab('baseline resistance') +
#   ylab('number of resistant drugs in treatment regimen') + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "lightgrey"),
#         panel.border = element_blank())
# 
# data %>% 
#   group_by(tx_resistant_drugs) %>% 
#   summarise(n = n()) %>% 
#   mutate(per = round((n / sum(n))*100,0))
# 
# # 7. create binaries for baseline resistance to each drug----
# subdata <- data %>% 
#   select(strain_id:Pas, num_resistant_drugs, effective_drugs, num_effective_drugs, tx_effective_drugs, tx_resistant_drugs)
# 
# drugs_short <- c('H', 'R', 'S', 'E', 'Z', 'Fq', 'Cm', 'Km', 'Am', 'Cs', 'Amx', 'Dld', 'Bdq', 'Ipm', 'Lzd', 'Cfz', 'Clr', 'Eto', 'Trd', 'Pas')
# 
# for (drug in drugs_short) {
#   val <- tb_drugs$get(drug)
#   subdata[,paste0('resistant_', drug)] <- case_when(grepl(val, subdata$resistant_drug, fixed = T) ~ 1,
#                                                      T ~ 0)
# }
# 
# subdata %>% 
#   select(H:Pas) %>% 
#   colSums()
# 
# for (drug in drugs_short) {
#   print(drug)
#   print(table(subdata[,drug], subdata[,paste0('resistant_', drug)]))
# }

# 8. create binary that indicates a mutation on ANY locus associated with drug resistance----
data <- data %>% 
  mutate(any_mutation = case_when(MTB000019 == 1 |
                                  MTB000020 == 1 |
                                  Rv0005 == 1 |
                                  Rv0006 == 1 |
                                  Rv0407 == 1 |
                                  Rv0486 == 1 |
                                  Rv0529 == 1 |
                                  Rv0667 == 1 |
                                  Rv0668 == 1 |
                                  Rv0676c == 1 |
                                  Rv0677c == 1 |
                                  Rv0678 == 1 |
                                  Rv0682 == 1 |
                                  Rv0701 == 1 |
                                  Rv0957 == 1 |
                                  Rv1173 == 1 |
                                  Rv1258c == 1 |
                                  Rv1267c == 1 |
                                  Rv1305 == 1 |
                                  Rv1484 == 1 |
                                  Rv1854c == 1 |
                                  Rv1908c == 1 |
                                  Rv1918c == 1 |
                                  Rv1979c == 1 |
                                  Rv2043c == 1 |
                                  Rv2416c == 1 |
                                  Rv2428 == 1 |
                                  Rv2535c == 1 |
                                  Rv2752c == 1 |
                                  Rv2780 == 1 |
                                  Rv2983 == 1 |
                                  Rv3083 == 1 |
                                  Rv3106 == 1 |
                                  Rv3197A == 1 |
                                  Rv3236c == 1 |
                                  Rv3261 == 1 |
                                  Rv3262 == 1 |
                                  Rv3423c == 1 |
                                  Rv3457c == 1 |
                                  Rv3547 == 1 |
                                  Rv3596c == 1 |
                                  Rv3601c == 1 |
                                  Rv3792 == 1 |
                                  Rv3793 == 1 |
                                  Rv3794 == 1 |
                                  Rv3795 == 1 |
                                  Rv3805c == 1 |
                                  Rv3806c == 1 |
                                  Rv3854c == 1 |
                                  Rv3855 == 1 |
                                  Rv3862c == 1 |
                                  Rv3919c == 1 ~ 1,
                                  T ~ 0))


# 9. save dataset with num of effective drugs----
saveRDS(data, file = "datasets/tx_hetmut_tbprof_neffectdrugs.Rds")
