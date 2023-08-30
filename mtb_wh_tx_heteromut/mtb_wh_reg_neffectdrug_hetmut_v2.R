#-------------------------------------------------------------------------#
# MSc Health Data Science 2023                                            #
# Dataset: Mtb Within-host genetic & tx data                              #
# Regressions: Predict hetero mutation based on number of effective drugs # 
#-------------------------------------------------------------------------#

library(tidyverse)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_hetmut_tbprof_neffectdrugs.Rds')

# 2. create count of total number of drugs treated with----
data <- data %>% 
  mutate(tx_total_drugs = rowSums((data[,14:34])))

# 3. simple model----
model_simple <- glm(any_mutation ~ tx_effective_drugs + tx_total_drugs, data = data, family = binomial(link = 'logit'))
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

# extract coefficient on total number of drugs
cat3 <- rownames(coef(summary(model_simple)))[3]
coef3 <- coef(summary(model_simple))[3,1]
pval3 <- coef(summary(model_simple))[3,4]
lower_ci3 <- confint(model_simple)[3,1]
upper_ci3 <- confint(model_simple)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

var_names <- c('predictor', 'log_estimate', 'log_lower_ci', 'log_upper_ci', 'pval')
results_simple <- rbind(results_simple, row1, row2, row3)

colnames(results_simple) <- var_names

results_simple <- results_simple %>% 
  mutate(estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(predictor, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(results_simple, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_reg_results_v2.csv')

# 4. adding in main_lineage----
# recode main_lineage into three categories: lineage 2, lineage 4, other (1, 3, 5)
data <- data %>% 
  mutate(main_lineage2 = case_when(main_lineage == 'lineage2' ~ 1,
                                   T ~ 0),
         main_lineage4 = case_when(main_lineage == 'lineage4' ~ 1,
                                   T ~ 0),
         main_lineageother = case_when(main_lineage != 'lineage2' & main_lineage != 'lineage4' ~ 1,
                                       T ~ 0))

model_lineage <- glm(any_mutation ~ tx_effective_drugs + tx_total_drugs + main_lineage2 + main_lineage4, data = data, family = binomial(link = 'logit'))
summary(model_lineage)

results_lineage <- data.frame()

# extract intercept
cat1 <- rownames(coef(summary(model_lineage)))[1]
coef1 <- coef(summary(model_lineage))[1,1]
pval1 <- coef(summary(model_lineage))[1,4]
lower_ci1 <- confint(model_lineage)[1,1]
upper_ci1 <- confint(model_lineage)[1,2]
row1 <- c(cat1, coef1, lower_ci1, upper_ci1, pval1)

# extract coefficient on number of effective drugs
cat2 <- rownames(coef(summary(model_lineage)))[2]
coef2 <- coef(summary(model_lineage))[2,1]
pval2 <- coef(summary(model_lineage))[2,4]
lower_ci2 <- confint(model_lineage)[2,1]
upper_ci2 <- confint(model_lineage)[2,2]
row2 <- c(cat2, coef2, lower_ci2, upper_ci2, pval2)

# extract coefficients total number of drugs
cat3 <- rownames(coef(summary(model_lineage)))[3] 
coef3 <- coef(summary(model_lineage))[3,1]
pval3 <- coef(summary(model_lineage))[3,4]
lower_ci3 <- confint(model_lineage)[3,1]
upper_ci3 <- confint(model_lineage)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

# extract coefficients on main lineage
cat4 <- rownames(coef(summary(model_lineage)))[4] # lineage2
coef4 <- coef(summary(model_lineage))[4,1]
pval4 <- coef(summary(model_lineage))[4,4]
lower_ci4 <- confint(model_lineage)[4,1]
upper_ci4 <- confint(model_lineage)[4,2]
row4 <- c(cat4, coef4, lower_ci4, upper_ci4, pval4)

cat5 <- rownames(coef(summary(model_lineage)))[5] # lineage4
coef5 <- coef(summary(model_lineage))[5,1]
pval5 <- coef(summary(model_lineage))[5,4]
lower_ci5 <- confint(model_lineage)[5,1]
upper_ci5 <- confint(model_lineage)[5,2]
row5 <- c(cat5, coef5, lower_ci5, upper_ci5, pval5)

results_lineage <- rbind(results_lineage, row1, row2, row3, row4, row5)

colnames(results_lineage) <- var_names

results_lineage <- results_lineage %>% 
  mutate(estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(predictor, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(results_lineage, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_lineage_reg_results_v2.csv')

# adding in baseline resistance 
model_baseresist <- glm(any_mutation ~ tx_effective_drugs + tx_total_drugs + baseline_resistance_type, data = data, family = binomial(link = 'logit'))
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

# extract coefficients on total number of drugs 
cat3 <- rownames(coef(summary(model_baseresist)))[3] # 
coef3 <- coef(summary(model_baseresist))[3,1]
pval3 <- coef(summary(model_baseresist))[3,4]
lower_ci3 <- confint(model_baseresist)[3,1]
upper_ci3 <- confint(model_baseresist)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

# extract coefficients on baseline resistance
cat4 <- rownames(coef(summary(model_baseresist)))[4] #  mono-resistant
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

write.csv(results_baseresist, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_baseresist_reg_results_v2.csv')

# full model with both main lineage and baseline resistance
model_full <- glm(any_mutation ~ tx_effective_drugs + tx_total_drugs + main_lineage2 + main_lineage4 + baseline_resistance_type, data = data, family = binomial(link = 'logit'))
summary(model_full)

results_full <- data.frame()

# extract intercept
cat1 <- rownames(coef(summary(model_full)))[1]
coef1 <- coef(summary(model_full))[1,1]
pval1 <- coef(summary(model_full))[1,4]
lower_ci1 <- confint(model_full)[1,1]
upper_ci1 <- confint(model_full)[1,2]
row1 <- c(cat1, coef1, lower_ci1, upper_ci1, pval1)

# extract coefficient on number of effective drugs 
cat2 <- rownames(coef(summary(model_full)))[2]
coef2 <- coef(summary(model_full))[2,1]
pval2 <- coef(summary(model_full))[2,4]
lower_ci2 <- confint(model_full)[2,1]
upper_ci2 <- confint(model_full)[2,2]
row2 <- c(cat2, coef2, lower_ci2, upper_ci2, pval2)

# extract coefficients on total number of drugs
cat3 <- rownames(coef(summary(model_full)))[3]
coef3 <- coef(summary(model_full))[3,1]
pval3 <- coef(summary(model_full))[3,4]
lower_ci3 <- confint(model_full)[3,1]
upper_ci3 <- confint(model_full)[3,2]
row3 <- c(cat3, coef3, lower_ci3, upper_ci3, pval3)

# extract coefficients on main lineage
cat4 <- rownames(coef(summary(model_full)))[4] # lineage2
coef4 <- coef(summary(model_full))[4,1]
pval4 <- coef(summary(model_full))[4,4]
lower_ci4 <- confint(model_full)[4,1]
upper_ci4 <- confint(model_full)[4,2]
row4 <- c(cat4, coef4, lower_ci4, upper_ci4, pval4)

cat5 <- rownames(coef(summary(model_full)))[5] # lineage4
coef5 <- coef(summary(model_full))[5,1]
pval5 <- coef(summary(model_full))[5,4]
lower_ci5 <- confint(model_full)[5,1]
upper_ci5 <- confint(model_full)[5,2]
row5 <- c(cat5, coef5, lower_ci5, upper_ci5, pval5)

# extract coefficients on baseline resistance 
cat6 <- rownames(coef(summary(model_full)))[6] # mono-resistant
coef6 <- coef(summary(model_full))[6,1]
pval6 <- coef(summary(model_full))[6,4]
lower_ci6 <- confint(model_full)[6,1]
upper_ci6 <- confint(model_full)[6,2]
row6 <- c(cat6, coef6, lower_ci6, upper_ci6, pval6)

cat7 <- rownames(coef(summary(model_full)))[7] # resistant
coef7 <- coef(summary(model_full))[7,1]
pval7 <- coef(summary(model_full))[7,4]
lower_ci7 <- confint(model_full)[7,1]
upper_ci7 <- confint(model_full)[7,2]
row7 <- c(cat7, coef7, lower_ci7, upper_ci7, pval7)

cat8 <- rownames(coef(summary(model_full)))[8] # mdr
coef8 <- coef(summary(model_full))[8,1]
pval8 <- coef(summary(model_full))[8,4]
lower_ci8 <- confint(model_full)[8,1]
upper_ci8 <- confint(model_full)[8,2]
row8 <- c(cat8, coef8, lower_ci8, upper_ci8, pval8)

cat9 <- rownames(coef(summary(model_full)))[9] # xdr
coef9 <- coef(summary(model_full))[9,1]
pval9 <- coef(summary(model_full))[9,4]
lower_ci9 <- confint(model_full)[9,1]
upper_ci9 <- confint(model_full)[9,2]
row9 <- c(cat9, coef9, lower_ci9, upper_ci9, pval9)

results_full <- rbind(results_full, row1, row2, row3, row4, row5, row6, row7, row8, row9)

colnames(results_full) <- var_names

results_full <- results_full %>% 
  mutate(estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(predictor, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(results_full, 'output/treatment_hetero_mutation/cont_neffectdrug_hetmut_allcovars_reg_results_v2.csv')
