#-------------------------------------------------------------------#
# MSc Health Data Science 2023                                      #
# Dataset: Mtb Within-host genetic & tx data                        #
# Regressions: Predict fixed mutation based on treatment   #
#-------------------------------------------------------------------#

library(tidyverse)
library(patchwork)

setwd("~/Desktop/LSHTM/Summer Project/analysis")

# 1. load data----
data <- readRDS('datasets/tx_hetmut_tbprof.Rds')

# 2. explore mutation data----
mutated_genes <- data %>% 
  select(MTB000001:Rv3924c)

num_mutated_genes <- mutated_genes %>% 
  colSums() %>% 
  data.frame() %>% 
  rename(num_strains = '.') %>% 
  arrange(desc(num_strains))

write.csv(num_mutated_genes, 'output/treatment_hetero_mutation/num_mutated_genes.csv')

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

write.csv(Rv0006_reg_results, 'output/treatment_hetero_mutation/Rv0006_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv0006_multi_reg1 <- glm(Rv0006 ~ H + R + E + Fq + Cm + Am + Cs + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg1)
# E has highest p-value 

Rv0006_multi_reg2 <- glm(Rv0006 ~ H + R + Fq + Cm + Am + Cs + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg2)
# H has highest p-value

Rv0006_multi_reg3 <- glm(Rv0006 ~ R + Fq + Cm + Am + Cs + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg3)
# Cs has highest p-value

Rv0006_multi_reg4 <- glm(Rv0006 ~ R + Fq + Cm + Am + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg4)
# Cm has highest p-value

Rv0006_multi_reg5 <- glm(Rv0006 ~ R + Fq + Am + Amx + Ipm + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg5)
# Eto has highest p-value

Rv0006_multi_reg6 <- glm(Rv0006 ~ R + Fq + Am + Amx + Ipm + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg6)
# Pas has highest p-value

Rv0006_multi_reg7 <- glm(Rv0006 ~ R + Fq + Am + Amx + Ipm + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg7)
# Amx has highest p-value

Rv0006_multi_reg8 <- glm(Rv0006 ~ R + Fq + Am + Ipm + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv0006_multi_reg8)
# R has highest p-value

Rv0006_multi_reg9 <- glm(Rv0006 ~ Fq + Am + Ipm + Cfz, data = data, family = binomial(link = 'logit'))
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

write.csv(Rv2043c_reg_results, 'output/treatment_hetero_mutation/Rv2043c_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv2043c_multi_reg1 <- glm(Rv2043c ~ H + R + Fq + Cm + Am + Cs + Cfz + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg1)
# Eto has highest p-value

Rv2043c_multi_reg2 <- glm(Rv2043c ~ H + R + Fq + Cm + Am + Cs + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg2)
# R has highest p-value

Rv2043c_multi_reg3 <- glm(Rv2043c ~ H + Fq + Cm + Am + Cs + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg3)
# H has highest p-value

Rv2043c_multi_reg4 <- glm(Rv2043c ~ Fq + Cm + Am + Cs + Cfz + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg4)
# Pas has highest p-value

Rv2043c_multi_reg5 <- glm(Rv2043c ~ Fq + Cm + Am + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg5)
# Cm has highest p-value

Rv2043c_multi_reg6 <- glm(Rv2043c ~ Fq + Am + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg6)
# Fq has highest p-value

Rv2043c_multi_reg7 <- glm(Rv2043c ~ Am + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg7)
# Am has highest p-value

Rv2043c_multi_reg8 <- glm(Rv2043c ~ Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv2043c_multi_reg8)
# final model

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

write.csv(Rv0667_reg_results, 'output/treatment_hetero_mutation/Rv0667_uni_results.csv')

# build multivariable model with drugs that were statistically significant
Rv0667_multi_reg1 <- glm(Rv0667 ~ S + Km, data = data, family = binomial(link = 'logit'))
summary(Rv0667_multi_reg1)
# final model

exp(coefficients(Rv0667_multi_reg1))
exp(confint(Rv0667_multi_reg1))
coef(summary(Rv0667_multi_reg1))[,4]

rm(Rv0667_multi_reg1)


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
write.csv(reg_results_wide, 'output/treatment_hetero_mutation/uni_results.csv')

# multi models on those with expected results
# 7. run backwards stepwise model on Rv3423c----
Rv3423c_multi_reg1 <- glm(Rv3423c ~ Cs + Amx + Ipm + Lzd + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg1)
# Pas has highest p-value

Rv3423c_multi_reg2 <- glm(Rv3423c ~ Cs + Amx + Ipm + Lzd + Eto, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg2)
# Amx has highest p-value

Rv3423c_multi_reg3 <- glm(Rv3423c ~ Cs + Ipm + Lzd + Eto, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg3)
# Lzd has highest p-value

Rv3423c_multi_reg4 <- glm(Rv3423c ~ Cs + Ipm + Eto, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg4)
# Cs has highest p-value

Rv3423c_multi_reg5 <- glm(Rv3423c ~ Ipm + Eto, data = data, family = binomial(link = 'logit'))
summary(Rv3423c_multi_reg5)
# final model

exp(coefficients(Rv3423c_multi_reg5))
exp(confint(Rv3423c_multi_reg5))
coef(summary(Rv3423c_multi_reg5))[,4]

rm(Rv3423c_multi_reg1,
   Rv3423c_multi_reg2,
   Rv3423c_multi_reg3,
   Rv3423c_multi_reg4,
   Rv3423c_multi_reg5)

# 8. run backwards stepwise model on Rv2780----
Rv2780_multi_reg1 <- glm(Rv2780 ~ E + Cm + Km + Cs + Amx + Ipm + Lzd + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg1)
# Lzd has highest p-value

Rv2780_multi_reg2 <- glm(Rv2780 ~ E + Cm + Km + Cs + Amx + Ipm + Eto + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg2)
# Eto has highest p-value

Rv2780_multi_reg3 <- glm(Rv2780 ~ E + Cm + Km + Cs + Amx + Ipm + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg3)
# Cs has highest p-value

Rv2780_multi_reg4 <- glm(Rv2780 ~ E + Cm + Km + Amx + Ipm + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg4)
# E has highest p-value

Rv2780_multi_reg5 <- glm(Rv2780 ~ Cm + Km + Amx + Ipm + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg5)
# Amx has highest p-value

Rv2780_multi_reg6 <- glm(Rv2780 ~ Cm + Km + Ipm + Pas, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg6)
# Pas has highest p-value

Rv2780_multi_reg7 <- glm(Rv2780 ~ Cm + Km + Ipm, data = data, family = binomial(link = 'logit'))
summary(Rv2780_multi_reg7)
# final model

exp(coefficients(Rv2780_multi_reg7))
exp(confint(Rv2780_multi_reg7))
coef(summary(Rv2780_multi_reg7))[,4]

rm(Rv2780_multi_reg1,
   Rv2780_multi_reg2,
   Rv2780_multi_reg3,
   Rv2780_multi_reg4,
   Rv2780_multi_reg5,
   Rv2780_multi_reg6,
   Rv2780_multi_reg7)

# 9. run backwards stepwise model on MTB000019----
MTB000019_multi_reg1 <- glm(MTB000019 ~ H + E + Z + Fq + Cm + Km + Cs + Amx + Bdq + Lzd + Cfz + Trd, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg1)
# Lzd has highest p-value

MTB000019_multi_reg2 <- glm(MTB000019 ~ H + E + Z + Fq + Cm + Km + Cs + Amx + Bdq + Cfz + Trd, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg2)
# Bdq has highest p-value

MTB000019_multi_reg3 <- glm(MTB000019 ~ H + E + Z + Fq + Cm + Km + Cs + Amx + Cfz + Trd, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg3)
# Trd has highest p-value

MTB000019_multi_reg4 <- glm(MTB000019 ~ H + E + Z + Fq + Cm + Km + Cs + Amx + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg4)
# H has highest p-value

MTB000019_multi_reg5 <- glm(MTB000019 ~ E + Z + Fq + Cm + Km + Cs + Amx + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg5)
# Z has highest p-value

MTB000019_multi_reg6 <- glm(MTB000019 ~ E + Fq + Cm + Km + Cs + Amx + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg6)
# Amx has highest p-value

MTB000019_multi_reg7 <- glm(MTB000019 ~ E + Fq + Cm + Km + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg7)
# Cm has highest p-value

MTB000019_multi_reg8 <- glm(MTB000019 ~ E + Fq + Km + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg8)
# Km has highest p-value

MTB000019_multi_reg8 <- glm(MTB000019 ~ E + Fq + Cs + Cfz, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg8)
# Cfz has highest p-value

MTB000019_multi_reg9 <- glm(MTB000019 ~ E + Fq + Cs, data = data, family = binomial(link = 'logit'))
summary(MTB000019_multi_reg9)
# final model

exp(coefficients(MTB000019_multi_reg9))
exp(confint(MTB000019_multi_reg9))
coef(summary(MTB000019_multi_reg9))[,4]

rm(MTB000019_multi_reg1,
   MTB000019_multi_reg2,
   MTB000019_multi_reg3,
   MTB000019_multi_reg4,
   MTB000019_multi_reg5,
   MTB000019_multi_reg6,
   MTB000019_multi_reg7,
   MTB000019_multi_reg8,
   MTB000019_multi_reg9)

# 10. run backwards stepwise model on Rv0678----
Rv0678_multi_reg1 <- glm(Rv0678 ~ S + Am + Dld + Bdq + Cfz, data = data, family = binomial(link = 'logit'))
summary(Rv0678_multi_reg1)
# Cfz has highest p-value

Rv0678_multi_reg2 <- glm(Rv0678 ~ S + Am + Dld + Bdq, data = data, family = binomial(link = 'logit'))
summary(Rv0678_multi_reg2)
# Bdq has highest p-value

Rv0678_multi_reg3 <- glm(Rv0678 ~ S + Am + Dld, data = data, family = binomial(link = 'logit'))
summary(Rv0678_multi_reg3)
# final model

exp(coefficients(Rv0678_multi_reg3))
exp(confint(Rv0678_multi_reg3))
coef(summary(Rv0678_multi_reg3))[,4]

rm(Rv0678_multi_reg1,
   Rv0678_multi_reg2,
   Rv0678_multi_reg3)
