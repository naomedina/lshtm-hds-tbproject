#-------------------------------------------------------------------------#
# MSc Health Data Science 2023                                            #
# Dataset: Mtb Within-host genetic & tx data                              #
# Merge in length of tx for those in TB Portals                           # 
#-------------------------------------------------------------------------#

library(tidyverse)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_hetmut_tbprof_neffectdrugs.Rds')
tbportals_tx <- read.csv('data/TBPortals/clinical_data/TB_Portals_Published_Clinical_Data_January_2023/TB_Portals_Regimens_January_2023.csv')

# 2. restrict data to samples in TB Portals----
data_tbportals <- data %>% 
  filter(study_id == 'tb_portals') %>% 
  select(strain_id:mutated_locus, num_resistant_drugs, effective_drugs, num_effective_drugs, tx_effective_drugs, tx_resistant_drugs, any_mutation)

# 3. merge in tx data from TB Portals and calculate total tx time----
data_tbportals <- data_tbportals %>% 
  merge(y = tbportals_tx, by = 'patient_id') %>% 
  group_by(patient_id) %>% 
  mutate(period_span = as.numeric(period_span),
         total_tx_time = case_when(n_distinct(condition_id) == 1 ~ mean(period_span),
                                   n_distinct(condition_id) > 1 ~ sum(period_span),
                                   T ~ NA),
         total_tx_time_yrs = total_tx_time / 365.25) %>% 
  select(patient_id:any_mutation, total_tx_time, total_tx_time_yrs) %>% 
  filter(!duplicated(strain_id)) %>% 
  ungroup()

# 4. summarise total tx time----
# summary(data_tbportals$total_tx_time)
# 
# data_tbportals %>% 
#   group_by(baseline_resistance_type) %>% 
#   summarise(min = min(total_tx_time, na.rm = T),
#             p25 = quantile(total_tx_time, 0.25, na.rm = T),
#             median = median(total_tx_time, na.rm = T),
#             mean = mean(total_tx_time, na.rm = T),
#             p75 = quantile(total_tx_time, 0.75, na.rm = T),
#             max = max(total_tx_time, na.rm = T))
# 
# hist(data_tbportals$total_tx_time,
#      col = 'cadetblue3',
#      main = 'total treatment time (days)',
#      xlab = '')
# 
# data_tbportals %>% ggplot(aes(x = baseline_resistance_type, y = total_tx_time)) +
#   geom_boxplot(fill = 'cadetblue3') +
#   xlab('baseline resistance') +
#   ylab('total treatment time (days)') + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "lightgrey"),
#         panel.border = element_blank())
# 
# data_tbportals %>% ggplot(aes(x = baseline_resistance_type, y = total_tx_time_yrs)) +
#   geom_boxplot(fill = 'cadetblue3') +
#   xlab('baseline resistance') +
#   ylab('total treatment time (days)') + 
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "lightgrey"),
#         panel.border = element_blank())
  
  
# 5. run analysis including total_tx_time on acquisition of any mutation----

# simplest model with only number of effective drugs + total tx time
model_simple <- glm(any_mutation ~ tx_effective_drugs + total_tx_time, data = data_tbportals, family = binomial(link = 'logit'))
summary(model_simple)

results_simple <- data.frame()

# extract intercept
cat1 <- rownames(coef(summary(model_simple)))[1]
coef1 <- coef(summary(model_simple))[1,1]
pval1 <- coef(summary(model_simple))[1,4]
lower_ci1 <- confint(model_simple)[1,1]
upper_ci1 <- confint(model_simple)[1,2]
row1 <- c(cat1, coef1, lower_ci1, upper_ci1, pval1)

# extract coefficient on number of effective drugs
cat2 <- rownames(coef(summary(model_simple)))[2]
coef2 <- coef(summary(model_simple))[2,1]
pval2 <- coef(summary(model_simple))[2,4]
lower_ci2 <- confint(model_simple)[2,1]
upper_ci2 <- confint(model_simple)[2,2]
row2 <- c(cat2, coef2, lower_ci2, upper_ci2, pval2)

# extract coefficient on length of tx
cat3 <- rownames(coef(summary(model_simple)))[3]
coef3 <- coef(summary(model_simple))[3,1]
pval3 <- coef(summary(model_simple))[3,4]
lower_ci3 <- confint(model_simple)[3,1]
upper_ci3 <- confint(model_simple)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

results_simple <- rbind(results_simple, row1, row2, row3)

var_names <- c('predictor', 'log_estimate', 'log_lower_ci', 'log_upper_ci', 'pval')
colnames(results_simple) <- var_names

results_simple <- results_simple %>% 
  mutate(estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(predictor, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(results_simple, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_txtime_reg_results.csv')

# adding in baseline resistance
model_baseresist <- glm(any_mutation ~ tx_effective_drugs + total_tx_time + baseline_resistance_type, data = data_tbportals, family = binomial(link = 'logit'))
summary(model_baseresist)

results_baseresist <- data.frame()

# extract intercept
cat1 <- rownames(coef(summary(model_baseresist)))[1]
coef1 <- coef(summary(model_baseresist))[1,1]
pval1 <- coef(summary(model_baseresist))[1,4]
lower_ci1 <- confint(model_baseresist)[1,1]
upper_ci1 <- confint(model_baseresist)[1,2]
row1 <- c(cat1, coef1, lower_ci1, upper_ci1, pval1)

# extract coefficient on number of effective drugs 
cat2 <- rownames(coef(summary(model_baseresist)))[2]
coef2 <- coef(summary(model_baseresist))[2,1]
pval2 <- coef(summary(model_baseresist))[2,4]
lower_ci2 <- confint(model_baseresist)[2,1]
upper_ci2 <- confint(model_baseresist)[2,2]
row2 <- c(cat2, coef2, lower_ci2, upper_ci2, pval2)

# extract coefficients on total tx time
cat3 <- rownames(coef(summary(model_baseresist)))[3] # lineage2
coef3 <- coef(summary(model_baseresist))[3,1]
pval3 <- coef(summary(model_baseresist))[3,4]
lower_ci3 <- confint(model_baseresist)[3,1]
upper_ci3 <- confint(model_baseresist)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

# extract coefficients on baseline resistance 
cat4 <- rownames(coef(summary(model_baseresist)))[4] # mono-resistant
coef4 <- coef(summary(model_baseresist))[4,1]
pval4 <- coef(summary(model_baseresist))[4,4]
lower_ci4 <- confint(model_baseresist)[4,1]
upper_ci4 <- confint(model_baseresist)[4,2]
row4 <- c(cat4, coef4, lower_ci4, upper_ci4, pval4)

cat5 <- rownames(coef(summary(model_baseresist)))[5] # resistant
coef5 <- coef(summary(model_baseresist))[5,1]
pval5 <- coef(summary(model_baseresist))[5,4]
lower_ci5 <- confint(model_baseresist)[5,1]
upper_ci5 <- confint(model_baseresist)[5,2]
row5 <- c(cat5, coef5, lower_ci5, upper_ci5, pval5)

cat6 <- rownames(coef(summary(model_baseresist)))[6] # mdr
coef6 <- coef(summary(model_baseresist))[6,1]
pval6 <- coef(summary(model_baseresist))[6,4]
lower_ci6 <- confint(model_baseresist)[6,1]
upper_ci6 <- confint(model_baseresist)[6,2]
row6 <- c(cat6, coef6, lower_ci6, upper_ci6, pval6)

cat7 <- rownames(coef(summary(model_baseresist)))[7] # xdr
coef7 <- coef(summary(model_baseresist))[7,1]
pval7 <- coef(summary(model_baseresist))[7,4]
lower_ci7 <- confint(model_baseresist)[7,1]
upper_ci7 <- confint(model_baseresist)[7,2]
row7 <- c(cat7, coef7, lower_ci7, upper_ci7, pval7)

results_baseresist <- rbind(results_baseresist, row1, row2, row3, row4, row5, row6, row7)

colnames(results_baseresist) <- var_names

results_baseresist <- results_baseresist %>% 
  mutate(estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(predictor, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(results_baseresist, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_baseresist_txtime_reg_results.csv')

