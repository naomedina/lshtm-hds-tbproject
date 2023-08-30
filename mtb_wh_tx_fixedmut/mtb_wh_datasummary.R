#-------------------------------------------------------------#
# MSc Health Data Science 2023                                #
# Dataset: Mtb Within-host genetic & tx data                  #
# Summary stats/plots                                         #
#-------------------------------------------------------------#

library(tidyverse)
library(patchwork)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_fixedmut_tbprof_neffectdrugs.Rds')

# 2. summary stats for table----
# full dataset: geography
data %>% 
  group_by(geography_recode) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0))
# stratified by baseline resistance type: geography
data %>% 
  group_by(baseline_resistance_type, geography_recode) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0)) %>% 
  print(n = Inf)
# stratified by any mutation
data %>% 
  group_by(any_mutation, geography_recode) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0)) %>% 
  print(n = Inf)

# full dataset: main lineage
data %>% 
  group_by(main_lineage) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0))
# stratified by baseline resistance type: main lineage
data %>% 
  group_by(baseline_resistance_type, main_lineage) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0)) %>% 
  print(n = Inf)
# stratified by any mutation
data %>% 
  group_by(any_mutation, main_lineage) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0)) %>% 
  print(n = Inf) 

# full dataset: baseline resistance
data %>% 
  group_by(baseline_resistance_type) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0))
# stratified by any mutation
data %>% 
  group_by(any_mutation, baseline_resistance_type) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100, 0)) %>% 
  print(n = Inf) 

# full dataset: treatment
data %>%
  select(H:Amx, Dld:Pas, baseline_resistance_type) %>%
  pivot_longer(
    cols = (H:Pas),
    names_to = 'drug',
    values_to = 'val'
  ) %>%
  filter(val==1) %>%
  group_by(drug) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 550,
         per = round(prop*100,0)) %>%
  print(n = Inf)

# stratified by any mutation: treatment
data %>%
  select(H:Amx, Dld:Pas, any_mutation) %>%
  pivot_longer(
    cols = (H:Pas),
    names_to = 'drug',
    values_to = 'val'
  ) %>% 
  arrange(any_mutation) %>% 
  filter(val == 1, any_mutation == 0) %>%
  group_by(any_mutation, drug) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 467,
         per = round(prop*100,0)) %>%
  print(n = Inf)

data %>%
  select(H:Amx, Dld:Pas, any_mutation) %>%
  pivot_longer(
    cols = (H:Pas),
    names_to = 'drug',
    values_to = 'val'
  ) %>% 
  arrange(any_mutation) %>% 
  filter(val == 1, any_mutation == 1) %>%
  group_by(any_mutation, drug) %>%
  summarise(n = n()) %>%
  mutate(prop = n / 83,
         per = round(prop*100,0)) %>%
  print(n = Inf)

# number of effective drugs
data %>% 
  group_by(baseline_resistance_type, tx_effective_drugs) %>% 
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n),
         per = round(prop*100,0)) %>% 
  print(n = Inf)

data %>% 
  group_by(any_mutation) %>% 
  summarise(min = min(tx_effective_drugs),
            p25 = quantile(tx_effective_drugs, .25),
            median = median(tx_effective_drugs),
            mean = mean(tx_effective_drugs),
            p75 = quantile(tx_effective_drugs, .75),
            max = max(tx_effective_drugs))

data %>% 
  summarise(min = min(tx_effective_drugs),
            p25 = quantile(tx_effective_drugs, .25),
            median = median(tx_effective_drugs),
            mean = mean(tx_effective_drugs),
            p75 = quantile(tx_effective_drugs, .75),
            max = max(tx_effective_drugs),
            )

# 11. plots----
# geography
geo_plot <- data %>% ggplot(aes(x = geography_recode)) + 
  geom_bar(fill = 'cadetblue3') + 
  xlab('') + 
  ylim(0,300) +
  theme_bw() + 
  ggtitle('2a. Geography') +
  theme(axis.line = element_line(colour = 'lightgrey'), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")) 

# main lineage
lineage_plot <- data %>% ggplot(aes(x = main_lineage)) + 
  geom_bar(fill = 'darkolivegreen3') + 
  xlab('') + 
  ylim(0, 300) +
  theme_bw() + 
  ggtitle('2b. Main lineage') + 
  theme(axis.line = element_line(colour = 'lightgrey'), 
        panel.border = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold")) 


# baseline resistance profile
# recode baseline resistance type for purpose of the plot

data$baseline_resistance_type <- recode(data$baseline_resistance_type, 
                                        'pan-susceptible' = 'pan-susceptible', 
                                        'mono-resisant' = 'mono-resistant',
                                        'resistant' = 'other resistant*',
                                        'MDR' = 'MDR',
                                        'XDR' = 'XDR',)
data$baseline_resistance_type <- factor(data$baseline_resistance_type,
                                        levels = c('pan-susceptible', 'mono-resistant', 'MDR', 'XDR', 'other resistant*'))

resist_plot <- data %>% ggplot(aes(x = baseline_resistance_type)) +
  geom_bar(fill = 'plum4') + 
  xlab('') + 
  ylim(0, 300) +
  theme_bw() + 
  ggtitle('2c. Baseline drug resistance') +
  theme(axis.line = element_line(colour = 'lightgrey'), 
        panel.border = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(size = 15, face = "bold"))


# plot treatment data
txplot <- data %>% 
  select(H:Pas) %>% 
  pivot_longer(
    cols = (H:Pas),
    names_to = 'drug',
    values_to = 'val'
  ) %>% 
  mutate(drug = factor(drug, 
                       levels = c('H', 'R', 'S', 'E', 'Z', 'Fq', 'Cm', 'Km', 'Am', 'Mfx', 'Cs', 'Amx', 'Dld', 'Bdq', 'Ipm', 'Lzd', 'Cfz', 'Clr', 'Eto', 'Trd', 'Pas'), 
                       labels = c('Isoniazid', 'Rifampicin', 'Streptomycin', 'Ethambutol', 'Pyrazinamide', 'Fluoroquinolones', 'Capreomycin', 'Kanamycin', 'Amikacin', 'Moxifloxacin', 'Cycloserine', 'Amoxicillin-clavulanate', 'Delamanid', 'Bedaquiline', 'Imipenem-cilastatin', 'Linezolid', 'Clofazimine', 'Clarithromycin', 'Pro/Ethionamide', 'Terizidone', 'Para-aminosalicylic acid'))) %>% 
  filter(val == 1, drug != 'Mb') %>% 
  ggplot(aes(x = drug)) + 
  geom_bar(fill = 'lightblue') +
  xlab('') + 
  theme_bw() + 
  ggtitle('2d. Treatment') +
  theme(axis.line = element_line(colour = 'lightgrey'), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, face = "bold")) 

plot <- geo_plot / lineage_plot / resist_plot / txplot
#plot <- plot + plot_annotation(title = 'Figure X. Dataset Characteristics')

#ggsave(filename = 'raw_plots/characteristics.jpeg', plot = plot, width=15, height=9)
#ggsave(filename = 'raw_plots/tx_summary.jpeg', plot = txplot, width=15, height=9)

ggsave(filename = 'raw_plots/all_chars.jpeg', plot = plot, width = 15, height = 18)


# plot of number of effective drugs by baseline resistance type
n_drugs <- data %>%
  group_by(baseline_resistance_type, tx_effective_drugs) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  data.frame() 
  
n_drugs_app <- data.frame(baseline_resistance_type = rep(c('pan-susceptible', 'mono-resistant', 'MDR', 'XDR', 'other resistant*'), times = 12),
                          tx_effective_drugs = rep(seq(0, 11), times = 5),
                          n = 0) 

n_drugs <- n_drugs %>% 
  merge(y = n_drugs_app, by = c('baseline_resistance_type', 'tx_effective_drugs'), all.y = T) %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  select(baseline_resistance_type:n.x) %>% 
  rename(n = n.x)

n_drugs_plot <- n_drugs %>%   
  ggplot(aes(x = tx_effective_drugs, y = n, fill = baseline_resistance_type)) +
  geom_bar(position = "dodge", stat = "identity") + 
  xlab('number of effective drugs used during treatment') + 
  ylab('number of strains') + 
  theme_bw() + 
  theme(axis.line = element_line(colour = 'lightgrey'), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_x_continuous(breaks = seq(0,11,1)) + 
  scale_fill_manual(name = 'Baseline drug resistance', 
                     breaks = c('pan-susceptible', 'mono-resistant', 'MDR','XDR', 'other resistant*'),
                     values = c('pan-susceptible'="cornflowerblue", 'mono-resistant'="cyan3", 'MDR'="dodgerblue3", 'XDR'="lightsteelblue3", 'other resistant*'="cadetblue2"))

ggsave(filename = 'raw_plots/num_effective_drugs.jpeg', plot = n_drugs_plot, width = 15, height = 9)
