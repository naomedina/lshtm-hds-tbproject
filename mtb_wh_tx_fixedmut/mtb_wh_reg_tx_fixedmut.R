#-------------------------------------------------------------------#
# MSc Health Data Science 2023                                      #
# Dataset: Mtb Within-host genetic & tx data                        #
# Regressions: Predict fixed mutation based on treatment   #
#-------------------------------------------------------------------#

library(tidyverse)
library(patchwork)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_fixedmut_tbprof.Rds')

# 2. explore mutation data----
mutated_genes <- data %>% 
  select(MTB000001:Rv3924c)

num_mutated_genes <- mutated_genes %>% 
  colSums() %>% 
  data.frame() %>% 
  rename(num_strains = '.') %>% 
  arrange(desc(num_strains))

write.csv(num_mutated_genes, 'output/treatment_fixed_mutation/num_mutated_genes.csv')

# 3. run backwards stepwise model on Rv0006 (gyrA: associated with resistance to fluoroquinolones)----
# univariate regs
drugs <- c('H', 'R', 'S', 'E', 'Z', 'Fq', 'Cm', 'Km', 'Am', 'Cs', 'Amx', 'Dld', 'Bdq', 'Ipm', 'Lzd', 'Cfz', 'Clr', 'Eto', 'Trd', 'Pas')

Rv0006_reg_results <- data.frame()

for (drug in drugs) {
  model <- glm(as.formula(paste0('Rv0006', '~', drug)), data = data, family = binomial(link = 'logit'))
  coef <- coef(summary(model))[2,1]
  pval <- coef(summary(model))[2,4]
  lower_ci <- confint(model)[2,1]
  upper_ci <- confint(model)[2,2]
  row <- c(drug, 'Rv0006', coef, lower_ci, upper_ci, pval)
  Rv0006_reg_results <- rbind(Rv0006_reg_results, row)
}

var_names <- c('drug', 'gene', 'log_estimate', 'log_lower_ci', 'log_upper_ci', 'pval')
colnames(Rv0006_reg_results) <- var_names

Rv0006_reg_results <- Rv0006_reg_results %>% 
  mutate(gene = factor(gene),
         estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(drug, gene, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(Rv0006_reg_results, 'output/treatment_fixed_mutation/Rv0006_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv0006_multi_reg1 <- glm(Rv0006 ~ H + R + E + Z + Fq + Cm + Am + Cs + Amx + Ipm + Lzd + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg1)
# H has highest p-value 

Rv0006_multi_reg2 <- glm(Rv0006 ~ R + E + Z + Fq + Cm + Am + Cs + Amx + Ipm + Lzd + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg2)
# E has highest p-value

Rv0006_multi_reg3 <- glm(Rv0006 ~ R + Z + Fq + Cm + Am + Cs + Amx + Ipm + Lzd + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg3)
# Cs has highest p-value

Rv0006_multi_reg4 <- glm(Rv0006 ~ R + Z + Fq + Cm + Am + Amx + Ipm + Lzd + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg4)
# Cm has highest p-value

Rv0006_multi_reg5 <- glm(Rv0006 ~ R + Z + Fq + Am + Amx + Ipm + Lzd + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg5)
# Lzd has highest p-value

Rv0006_multi_reg6 <- glm(Rv0006 ~ R + Z + Fq + Am + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg6)
# Eto has highest p-value

Rv0006_multi_reg7 <- glm(Rv0006 ~ R + Z + Fq + Am + Amx + Ipm + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg7)
# Amx has highest p-value

Rv0006_multi_reg8 <- glm(Rv0006 ~ R + Z + Fq + Am + Ipm + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg8)
# Pas has highest p-value

Rv0006_multi_reg9 <- glm(Rv0006 ~ R + Z + Fq + Am + Ipm + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg9)
# final model

exp(coefficients(Rv0006_multi_reg9))
exp(confint(Rv0006_multi_reg9))
coef(summary(Rv0006_multi_reg9))[,4]

rm(Rv0006_multi_reg1, 
   Rv0006_multi_reg2, 
   Rv0006_multi_reg3, 
   Rv0006_multi_reg4, 
   Rv0006_multi_reg5, 
   Rv0006_multi_reg6, 
   Rv0006_multi_reg7, 
   Rv0006_multi_reg8, 
   Rv0006_multi_reg9)

# plot results of univariate models
# Rv0006_plot_center <- Rv0006_reg_results %>% 
#   ggplot(aes(y = fct_rev(drug))) + 
#   theme_classic() +
#   geom_point(aes(x = log_estimate), shape = 15, size = 3) +
#   geom_linerange(aes(xmin = log_lower_ci, xmax = log_upper_ci)) +
#   geom_vline(xintercept = 0, linetype = 'dashed') +
#   labs(x = 'Log-Odds Ratio', y = '') +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) +
#   coord_cartesian(ylim=c(1,20))
# 
# Rv0006_plot_data <- Rv0006_reg_results %>% 
#   mutate(across(c(estimate, lower_ci, upper_ci),
#                 ~ format(round(.x, digits = 2), nsmall = 2)),
#          estimate_lab = paste0(estimate, ' (', lower_ci, '-', upper_ci, ')')) %>% 
#   mutate(pval = case_when(pval < .001 ~ '<0.001',
#                           round(pval, 2) == 0.05 ~ as.character(round(pval, 3)),
#                           pval < .01 ~ str_pad(as.character(round(pval, 3)),
#                                                width = 4,
#                                                pad = '0',
#                                                side = 'right'),
#                           T ~ str_pad(as.character(round(pval, 2)),
#                                       width = 4,
#                                       pad = '0',
#                                       side = 'right'))) %>% 
#   bind_rows(data.frame(drug = 'Drug', 
#                        estimate_lab = 'Odds Ratio (95% CI)',
#                        lower_ci = '',
#                        upper_ci = '',
#                        pval = 'p-value')) %>% 
#   mutate(model = fct_rev(fct_relevel(drug, 'Drug')))
# 
# 
# Rv0006_plot_left <- Rv0006_plot_data %>% 
#   ggplot(aes(y = model)) +
#   geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = estimate_lab),
#     hjust = 0,
#     fontface = ifelse(Rv0006_plot_data$estimate_lab == "Odds Ratio (95% CI)", "bold", "plain")) +
#   theme_void() +
#   coord_cartesian(xlim = c(0,4))
# 
# Rv0006_plot_right <- Rv0006_plot_data %>% 
#   ggplot() + 
#   geom_text(
#     aes(x = 0, y = model, label = pval),
#     hjust = 0,
#     fontface = ifelse(Rv0006_plot_data$pval == 'p-value', 'bold', 'plain')) +
#   theme_void()
# 
# layout <- c(
#   area(t = 0, l = 0, b = 30, r = 3),
#   area(t = 2, l = 4, b = 30, r = 9),
#   area(t = 0, l = 9, b = 30, r = 11)
# )
# 
# Rv0006_plot_left + Rv0006_plot_center + Rv0006_plot_right + plot_layout(design = layout)

# 4. run backwards stepwise model on Rv2043c (pncA: associated with resistance to pyrazinamide)----
Rv2043c_reg_results <- data.frame()

for (drug in drugs) {
  model <- glm(as.formula(paste0('Rv2043c', '~', drug)), data = data, family = binomial(link = 'logit'))
  coef <- coef(summary(model))[2,1]
  pval <- coef(summary(model))[2,4]
  lower_ci <- confint(model)[2,1]
  upper_ci <- confint(model)[2,2]
  row <- c(drug, 'Rv2043c', coef, lower_ci, upper_ci, pval)
  Rv2043c_reg_results <- rbind(Rv2043c_reg_results, row)
}

colnames(Rv2043c_reg_results) <- var_names

Rv2043c_reg_results <- Rv2043c_reg_results %>% 
  mutate(gene = factor(gene),
         estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(drug, gene, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(Rv2043c_reg_results, 'output/treatment_fixed_mutation/Rv2043c_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv2043c_multi_reg1 <- glm(Rv2043c ~ H + R + E + Cm + Am + Cs + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg1)
# R has highest p-value

Rv2043c_multi_reg2 <- glm(Rv2043c ~ H + E + Cm + Am + Cs + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg2)
# Eto has highest p-value

Rv2043c_multi_reg3 <- glm(Rv2043c ~ H + E + Cm + Am + Cs + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg3)
# Cm has highest p-value

Rv2043c_multi_reg4 <- glm(Rv2043c ~ H + E + Am + Cs + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg4)
# H has highest p-value

Rv2043c_multi_reg5 <- glm(Rv2043c ~ E + Am + Cs + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg5)
# E has highest p-value

Rv2043c_multi_reg6 <- glm(Rv2043c ~ Am + Cs + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg6)
# Am has highest p-value

Rv2043c_multi_reg7 <- glm(Rv2043c ~ Cs + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg7)
# Cs has highest p-value

Rv2043c_multi_reg8 <- glm(Rv2043c ~ Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg8)
# finalmodel

exp(coefficients(Rv2043c_multi_reg8))
exp(confint(Rv2043c_multi_reg8))
coef(summary(Rv2043c_multi_reg8))[,4]

rm(Rv2043c_multi_reg1, 
   Rv2043c_multi_reg2, 
   Rv2043c_multi_reg3, 
   Rv2043c_multi_reg4, 
   Rv2043c_multi_reg5, 
   Rv2043c_multi_reg6, 
   Rv2043c_multi_reg7, 
   Rv2043c_multi_reg8)

# plot results of univariate models
# Rv2043c_plot_center <- Rv2043c_reg_results %>% 
#   ggplot(aes(y = fct_rev(drug))) + 
#   theme_classic() +
#   geom_point(aes(x = log_estimate), shape = 15, size = 3) +
#   geom_linerange(aes(xmin = log_lower_ci, xmax = log_upper_ci)) +
#   geom_vline(xintercept = 0, linetype = 'dashed') +
#   labs(x = 'Log-Odds Ratio', y = '') +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) +
#   coord_cartesian(ylim=c(1,20))
# 
# Rv2043c_plot_data <- Rv2043c_reg_results %>% 
#   mutate(across(c(estimate, lower_ci, upper_ci),
#                 ~ format(round(.x, digits = 2), nsmall = 2)),
#          estimate_lab = paste0(estimate, ' (', lower_ci, '-', upper_ci, ')')) %>% 
#   mutate(pval = case_when(pval < .001 ~ '<0.001',
#                           round(pval, 2) == 0.05 ~ as.character(round(pval, 3)),
#                           pval < .01 ~ str_pad(as.character(round(pval, 3)),
#                                                width = 4,
#                                                pad = '0',
#                                                side = 'right'),
#                           T ~ str_pad(as.character(round(pval, 2)),
#                                       width = 4,
#                                       pad = '0',
#                                       side = 'right'))) %>% 
#   bind_rows(data.frame(drug = 'Drug', 
#                        estimate_lab = 'Odds Ratio (95% CI)',
#                        lower_ci = '',
#                        upper_ci = '',
#                        pval = 'p-value')) %>% 
#   mutate(model = fct_rev(fct_relevel(drug, 'Drug')))
# 
# 
# Rv2043c_plot_left <- Rv2043c_plot_data %>% 
#   ggplot(aes(y = model)) +
#   geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = estimate_lab),
#     hjust = 0,
#     fontface = ifelse(Rv2043c_plot_data$estimate_lab == "Odds Ratio (95% CI)", "bold", "plain")) +
#   theme_void() +
#   coord_cartesian(xlim = c(0,4))
# 
# Rv2043c_plot_right <- Rv2043c_plot_data %>% 
#   ggplot() + 
#   geom_text(
#     aes(x = 0, y = model, label = pval),
#     hjust = 0,
#     fontface = ifelse(Rv2043c_plot_data$pval == 'p-value', 'bold', 'plain')) +
#   theme_void()
# 
# layout <- c(
#   area(t = 2.5, l = 0, b = 30, r = 3),
#   area(t = 0, l = 4, b = 30, r = 9),
#   area(t = 2.5, l = 9, b = 30, r = 11)
# )
# 
# Rv2043c_plot_left + Rv2043c_plot_center + Rv2043c_plot_right + plot_layout(design = layout)

# 5. run backwards stepwise model on Rv0667 (rpoB: associated with resistance to rifampicin)----
Rv0667_reg_results <- data.frame()

for (drug in drugs) {
  model <- glm(as.formula(paste0('Rv0667', '~', drug)), data = data, family = binomial(link = 'logit'))
  coef <- coef(summary(model))[2,1]
  pval <- coef(summary(model))[2,4]
  lower_ci <- confint(model)[2,1]
  upper_ci <- confint(model)[2,2]
  row <- c(drug, 'Rv0667', coef, lower_ci, upper_ci, pval)
  Rv0667_reg_results <- rbind(Rv0667_reg_results, row)
}

colnames(Rv0667_reg_results) <- var_names

Rv0667_reg_results <- Rv0667_reg_results %>% 
  mutate(gene = factor(gene),
         estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval),
         sig_flag = case_when(pval <= 0.05 ~ 1,
                              T ~ 0)) %>% 
  select(drug, gene, estimate, lower_ci, upper_ci, pval, sig_flag)

write.csv(Rv0667_reg_results, 'output/treatment_fixed_mutation/Rv0667_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv0667_multi_reg1 <- glm(Rv0667 ~ Km + Eto, data = data, family = binomial(link = 'logit'))
summary(Rv0667_multi_reg1)
# Km has highest p-value

Rv0667_multi_reg2 <- glm(Rv0667 ~ Eto, data = data, family = binomial(link = 'logit'))
summary(Rv0667_multi_reg2)
# final model

exp(coefficients(Rv0667_multi_reg2))
exp(confint(Rv0667_multi_reg2))
coef(summary(Rv0667_multi_reg2))[,4]

rm(Rv0667_multi_reg1,
   Rv0667_multi_reg2) 

# plot results of univariate models
# Rv0667_plot_center <- Rv0667_reg_results %>% 
#   ggplot(aes(y = fct_rev(drug))) + 
#   theme_classic() +
#   geom_point(aes(x = log_estimate), shape = 15, size = 3) +
#   geom_linerange(aes(xmin = log_lower_ci, xmax = log_upper_ci)) +
#   geom_vline(xintercept = 0, linetype = 'dashed') +
#   labs(x = 'Log-Odds Ratio', y = '') +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y= element_blank(),
#         axis.text.y= element_blank(),
#         axis.title.y= element_blank()) +
#   coord_cartesian(ylim=c(1,20))
# 
# Rv0667_plot_data <- Rv0667_reg_results %>% 
#   mutate(across(c(estimate, lower_ci, upper_ci),
#                 ~ format(round(.x, digits = 2), nsmall = 2)),
#          estimate_lab = paste0(estimate, ' (', lower_ci, '-', upper_ci, ')')) %>% 
#   mutate(pval = case_when(pval < .001 ~ '<0.001',
#                           round(pval, 2) == 0.05 ~ as.character(round(pval, 3)),
#                           pval < .01 ~ str_pad(as.character(round(pval, 3)),
#                                                width = 4,
#                                                pad = '0',
#                                                side = 'right'),
#                           T ~ str_pad(as.character(round(pval, 2)),
#                                       width = 4,
#                                       pad = '0',
#                                       side = 'right'))) %>% 
#   bind_rows(data.frame(drug = 'Drug', 
#                        estimate_lab = 'Odds Ratio (95% CI)',
#                        lower_ci = '',
#                        upper_ci = '',
#                        pval = 'p-value')) %>% 
#   mutate(model = fct_rev(fct_relevel(drug, 'Drug')))
# 
# 
# Rv0667_plot_left <- Rv0667_plot_data %>% 
#   ggplot(aes(y = model)) +
#   geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold") +
#   geom_text(
#     aes(x = 1, label = estimate_lab),
#     hjust = 0,
#     fontface = ifelse(Rv0667_plot_data$estimate_lab == "Odds Ratio (95% CI)", "bold", "plain")) +
#   theme_void() +
#   coord_cartesian(xlim = c(0,4))
# 
# Rv0667_plot_right <- Rv0667_plot_data %>% 
#   ggplot() + 
#   geom_text(
#     aes(x = 0, y = model, label = pval),
#     hjust = 0,
#     fontface = ifelse(Rv0667_plot_data$pval == 'p-value', 'bold', 'plain')) +
#   theme_void()
# 
# layout <- c(
#   area(t = 0, l = 0, b = 30, r = 3),
#   area(t = 2.5, l = 4, b = 30, r = 9),
#   area(t = 0, l = 9, b = 30, r = 11)
# )
# 
# Rv0667_plot_left + Rv0667_plot_center + Rv0667_plot_right + plot_layout(design = layout)

# 6. Univariate regressions for remaining genes----
# loop for univariate regressions of remaining genes with relatively high number of strains
loci <- c('Rv3795', 'Rv1694', 'Rv3423c', 'Rv1908c', 'Rv2780', 'Rv1129c', 'Rv3854c', 'MTB000019', 'Rv0678', 'Rv0758', 'Rv0005', 'Rv1408', 'Rv2571c', 'Rv2764c', 'Rv3761c')

reg_results <- data.frame()

for (drug in drugs) {
  for (gene in loci) {
    model <- glm((as.formula(paste0(gene, '~', drug))), data = data, family = binomial(link = 'logit'))
    coef <- coef(summary(model))[2,1]
    pval <- coef(summary(model))[2,4]
    lower_ci <- confint(model)[2,1]
    upper_ci <- confint(model)[2,2]
    row <- c(drug, gene, coef, lower_ci, upper_ci, pval)
    reg_results <- rbind(reg_results, row)
  }
}

var_names <- c('drug', 'gene', 'log_estimate', 'log_lower_ci', 'log_upper_ci', 'pval')
colnames(reg_results) <- var_names

reg_results <- reg_results %>% 
  mutate(gene = factor(gene),
         estimate = exp(as.numeric(log_estimate)),
         lower_ci = exp(as.numeric(log_lower_ci)),
         upper_ci = exp(as.numeric(log_upper_ci)),
         pval = as.numeric(pval)) %>% 
  select(drug, gene, estimate, lower_ci, upper_ci, pval)

# reshape reg_results to be wide
reg_results_wide <- reg_results %>% 
  pivot_wider(
    id_cols = drug,
    names_from = gene,
    values_from = c(estimate, lower_ci, upper_ci, pval),
    names_vary = 'slowest',
    names_glue = '{gene}.{.value}'
  )

# export results
write.csv(reg_results_wide, 'output/treatment_fixed_mutation/uni_results.csv')

# note 31 JUL 2023: models from here have not been updated!
# 7. run backwards stepwise model on Rv3795----
# build multivariable model with drugs that were statistically significant
Rv3795_multi_reg1 <- glm(Rv3795 ~ Am + Km, data = data, family = binomial(link = 'logit'))
summary(Rv3795_multi_reg1)
# Eto has highest p-value

Rv3795_multi_reg2 <- glm(Rv3795 ~ Am + Km, data = data, family = binomial(link = 'logit'))
summary(Rv3795_multi_reg2)
# final model, all significant at 5% level

exp(coefficients(Rv3795_multi_reg2))
exp(confint(Rv3795_multi_reg2))
coef(summary(Rv3795_multi_reg2))[,4]

rm(Rv3795_multi_reg1,
   Rv3795_multi_reg2)

# 8. run backwards stepwise model on Rv3423c----
Rv3423c_multi_reg1 <- glm(Rv3423c ~ Ipm + Lzd, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg1)
# Ipm has highest p-value

Rv3423c_multi_reg2 <- glm(Rv3423c ~ Lzd, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg2)
# final model

exp(coefficients(Rv3423c_multi_reg2))
exp(confint(Rv3423c_multi_reg2))
coef(summary(Rv3423c_multi_reg2))[,4]

rm(Rv3423c_multi_reg1,
   Rv3423c_multi_reg2)

# 9. run backwards stepwise model on Rv2780----
Rv2780_multi_reg1 <- glm(Rv2780 ~ Cm + Mfx + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg1)
# Cm has highest p-value

Rv2780_multi_reg2 <- glm(Rv2780 ~ Mfx + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg2)
# Pas has highest p-value

Rv2780_multi_reg3 <- glm(Rv2780 ~ Mfx, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg3)
# final model

exp(coefficients(Rv2780_multi_reg3))
exp(confint(Rv2780_multi_reg3))
coef(summary(Rv2780_multi_reg3))[,4]

rm(Rv2780_multi_reg1,
   Rv2780_multi_reg2,
   Rv2780_multi_reg3)

# 10. run backwards stepwise model on Rv1129c----
Rv1129c_multi_reg1 <- glm(Rv1129c ~ Amx + Dld + Ipm + Lzd, data = data, family = binomial(link = 'logit'))
summary(Rv1129c_multi_reg1)
# Amx has highest p-value

Rv1129c_multi_reg2 <- glm(Rv1129c ~ Dld + Ipm + Lzd, data = data, family = binomial(link = 'logit'))
summary(Rv1129c_multi_reg2)
# Lzd has highest p-value

Rv1129c_multi_reg3 <- glm(Rv1129c ~ Dld + Ipm, data = data, family = binomial(link = 'logit'))
summary(Rv1129c_multi_reg3)
# Dld has highest p-value

Rv1129c_multi_reg4 <- glm(Rv1129c ~ Ipm, data = data, family = binomial(link = 'logit'))
summary(Rv1129c_multi_reg4)
# final model

exp(coefficients(Rv1129c_multi_reg4))
exp(confint(Rv1129c_multi_reg4))
coef(summary(Rv1129c_multi_reg4))[,4]

rm(Rv1129c_multi_reg1,
   Rv1129c_multi_reg2,
   Rv1129c_multi_reg3,
   Rv1129c_multi_reg4)

# 11. run backwards stepwise model on MTB000019----
MTB000019_multi_reg1 <- glm(MTB000019 ~ S + Am + Clr, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg1)
# S has highest p-value

MTB000019_multi_reg2 <- glm(MTB000019 ~ Am + Clr, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg2)
# final model

exp(coefficients(MTB000019_multi_reg2))
exp(confint(MTB000019_multi_reg2))
coef(summary(MTB000019_multi_reg2))[,4]

rm(MTB000019_multi_reg1,
   MTB000019_multi_reg2)

# 12. run backwards stepwise model on Rv0678----
Rv0678_multi_reg1 <- glm(Rv0678 ~ Dld + Bdq + Clr + Trd, data = data, family = binomial(link = 'logit'))
summary(Rv0678_multi_reg1)
# Trd has highest p-value

Rv0678_multi_reg2 <- glm(Rv0678 ~ Dld + Bdq + Clr, data = data, family = binomial(link = 'logit'))
summary(Rv0678_multi_reg2)
# final model

exp(coefficients(Rv0678_multi_reg2))
exp(confint(Rv0678_multi_reg2))
coef(summary(Rv0678_multi_reg2))[,4]

rm(Rv0678_multi_reg1,
   Rv0678_multi_reg2)

# 13. run backwards stepwise model on Rv0758----
Rv0758_multi_reg1 <- glm(Rv0758 ~ Amx + Ipm, data = data, family = binomial(link = 'logit'))
summary(Rv0758_multi_reg1)
# Amx has highest p-value

Rv0758_multi_reg2 <- glm(Rv0758 ~ Ipm, data = data, family = binomial(link = 'logit'))
summary(Rv0758_multi_reg2)
# final model

exp(coefficients(Rv0758_multi_reg2))
exp(confint(Rv0758_multi_reg2))
coef(summary(Rv0758_multi_reg2))[,4]

rm(Rv0758_multi_reg1,
   Rv0758_multi_reg2)
