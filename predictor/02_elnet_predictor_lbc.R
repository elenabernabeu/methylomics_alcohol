#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")

datadir <- "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/"
localdir <- "/Local_Data/methylation/GS_20k/Chromosomes/" # p17

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)


## Test in LBC
####################################################################

alcohol <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_10506_usualdrinkers.tsv")
alcohol_all <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_16717.tsv")
rownames(alcohol) <- alcohol$ID
rownames(alcohol_all) <- alcohol_all$ID

coefs_daniel <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/daniel_predictor.tsv", header = T)
coefs <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3w4_noadjustments_scaledmeth.tsv", header = T)
coefs_all <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_scaledmeth.tsv", header = T)
coefs_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)
coefs_all_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)

rownames(coefs_daniel) <- coefs_daniel$CpG
rownames(coefs) <- coefs$CpG_Site
rownames(coefs_all) <- coefs_all$CpG_Site
rownames(coefs_filt) <- coefs_filt$CpG_Site
rownames(coefs_all_filt) <- coefs_all_filt$CpG_Site

coefs_ni_daniel <- coefs_daniel[2:nrow(coefs_daniel),]
coefs_ni <- coefs[2:nrow(coefs),]
coefs_ni_all <- coefs_all[2:nrow(coefs_all),]
coefs_ni_filt <- coefs_filt[2:nrow(coefs_filt),]
coefs_ni_all_filt <- coefs_all_filt[2:nrow(coefs_all_filt),]

# Import LBC data
lbc_target <- read.table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T, row.names = 1)
lbc_mvals <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_mvals_3489.rds")
lbc_mvals <- t(lbc_mvals)
lbc_mvals <- m2beta(lbc_mvals)

# Prep pheno
lbc_target$alcunitsupw <- as.numeric(lbc_target$alcunitsupw)
lbc_target <- lbc_target[!is.na(lbc_target$alcunitsupw),]
lbc_target$alcunitsupw_log <- log(lbc_target$alcunitsupw + 1)

# Replace missing values with 0
lbc_mvals[is.na(lbc_mvals)] <- 0

# Divide LBC data into LBC21 and LBC36
lbc_target_36 <- lbc_target[lbc_target["cohort"] == "LBC36",] # 2796, waves 1, 2, 3, and 4
lbc_target_21 <- lbc_target[lbc_target["cohort"] == "LBC21",] # 692, waves 1, 3, and 4

# Keep just wave 1
lbc_target_21 <- lbc_target_21[lbc_target_21$WAVE == 1,] # 436
lbc_target_36 <- lbc_target_36[lbc_target_36$WAVE == 1,] # 895

# Match CpGs
coefs_ni_lbc <- coefs_ni[which(rownames(coefs_ni) %in% colnames(lbc_mvals)),] # 1035/1073
coefs_ni_lbc_all <- coefs_ni_all[which(rownames(coefs_ni_all) %in% colnames(lbc_mvals)),] # 1288/1327
coefs_ni_lbc_daniel <- coefs_ni_daniel[which(rownames(coefs_ni_daniel) %in% colnames(lbc_mvals)),] # 437/449
coefs_ni_lbc_filt <- coefs_ni_filt[which(rownames(coefs_ni_filt) %in% colnames(lbc_mvals)),] # 568/585 (315/325 Carreras)
coefs_ni_lbc_all_filt <- coefs_ni_all_filt[which(rownames(coefs_ni_all_filt) %in% colnames(lbc_mvals)),] # 647/659 (425/430)

lbc_mvals_daniel <- lbc_mvals[,rownames(coefs_ni_lbc_daniel)]
lbc_mvals_ud <- lbc_mvals[,rownames(coefs_ni_lbc)]
lbc_mvals_ud_filt <- lbc_mvals[,rownames(coefs_ni_lbc_filt)]
lbc_mvals_all <- lbc_mvals[,rownames(coefs_ni_lbc_all)]
lbc_mvals_all_filt <- lbc_mvals[,rownames(coefs_ni_lbc_all_filt)]

# Scale M-vals
#lbc_mvals_36 <- scale(lbc_mvals_36)
#lbc_mvals_21 <- scale(lbc_mvals_21)

# Subset
lbc_mvals_36 <- scale(lbc_mvals_ud[rownames(lbc_target_36),])
lbc_mvals_21 <- scale(lbc_mvals_ud[rownames(lbc_target_21),])
lbc_mvals_36_all <- scale(lbc_mvals_all[rownames(lbc_target_36),])
lbc_mvals_21_all <- scale(lbc_mvals_all[rownames(lbc_target_21),])
lbc_mvals_36_filt <- scale(lbc_mvals_ud_filt[rownames(lbc_target_36),])
lbc_mvals_21_filt <- scale(lbc_mvals_ud_filt[rownames(lbc_target_21),])
lbc_mvals_36_all_filt <- scale(lbc_mvals_all_filt[rownames(lbc_target_36),])
lbc_mvals_21_all_filt <- scale(lbc_mvals_all_filt[rownames(lbc_target_21),])
lbc_mvals_36_dan <- scale(lbc_mvals_daniel[rownames(lbc_target_36),])
lbc_mvals_21_dan <- scale(lbc_mvals_daniel[rownames(lbc_target_21),])

# Predictions
pred_36 <- lbc_mvals_36 %*% coefs_ni_lbc$Coefficient
pred_36 <- as.data.frame(pred_36)
names(pred_36) <- c("ac_pred")

pred_36_all <- lbc_mvals_36_all %*% coefs_ni_lbc_all$Coefficient
pred_36_all <- as.data.frame(pred_36_all)
names(pred_36_all) <- c("ac_pred")

pred_36_filt <- lbc_mvals_36_filt %*% coefs_ni_lbc_filt$Coefficient
pred_36_filt <- as.data.frame(pred_36_filt)
names(pred_36_filt) <- c("ac_pred")

pred_36_all_filt <- lbc_mvals_36_all_filt %*% coefs_ni_lbc_all_filt$Coefficient
pred_36_all_filt <- as.data.frame(pred_36_all_filt)
names(pred_36_all_filt) <- c("ac_pred")

pred_36_dan <- lbc_mvals_36_dan %*% coefs_ni_lbc_daniel$Beta
pred_36_dan <- as.data.frame(pred_36_dan)
names(pred_36_dan) <- c("ac_pred")

pred_21 <- lbc_mvals_21 %*% coefs_ni_lbc$Coefficient
pred_21 <- as.data.frame(pred_21)
names(pred_21) <- c("ac_pred")

pred_21_all <- lbc_mvals_21_all %*% coefs_ni_lbc_all$Coefficient
pred_21_all <- as.data.frame(pred_21_all)
names(pred_21_all) <- c("ac_pred")

pred_21_filt <- lbc_mvals_21_filt %*% coefs_ni_lbc_filt$Coefficient
pred_21_filt <- as.data.frame(pred_21_filt)
names(pred_21_filt) <- c("ac_pred")

pred_21_all_filt <- lbc_mvals_21_all_filt %*% coefs_ni_lbc_all_filt$Coefficient
pred_21_all_filt <- as.data.frame(pred_21_all_filt)
names(pred_21_all_filt) <- c("ac_pred")

pred_21_dan <- lbc_mvals_21_dan %*% coefs_ni_lbc_daniel$Beta
pred_21_dan <- as.data.frame(pred_21_dan)
names(pred_21_dan) <- c("ac_pred")

# Export
write.table(data.frame(basename = rownames(pred_36), pred_36), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21), pred_21), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_36_all), pred_36_all), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_all), pred_21_all), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_36_filt), pred_36_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_filt), pred_21_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_36_all_filt), pred_36_all_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_all_filt), pred_21_all_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_36_dan), pred_36_dan), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_danielpredictor_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_dan), pred_21_dan), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_danielpredictor_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)


## Import
####################################################################

pred_36 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", header = T, row.names = 1)
pred_21 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", header = T, row.names = 1)
pred_36_all <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", header = T, row.names = 1)
pred_21_all <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", header = T, row.names = 1)
pred_36_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T, row.names = 1)
pred_21_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T, row.names = 1)
pred_36_all_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T, row.names = 1)
pred_21_all_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T, row.names = 1)
pred_36_dan <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_danielpredictor_noadjustments.tsv", header = T, row.names = 1)
pred_21_dan <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_danielpredictor_noadjustments.tsv", header = T, row.names = 1)

lbc_target_36 <- lbc_target_36[rownames(pred_36),]
lbc_target_36$ac_pred <- pred_36$ac_pred
lbc_target_36$ac_pred_all <- pred_36_all$ac_pred
lbc_target_36$ac_pred_filt <- pred_36_filt$ac_pred
lbc_target_36$ac_pred_all_filt <- pred_36_all_filt$ac_pred
lbc_target_36$ac_pred_dan <- pred_36_dan$ac_pred

lbc_target_21 <- lbc_target_21[rownames(pred_21),]
lbc_target_21$ac_pred <- pred_21$ac_pred
lbc_target_21$ac_pred_all <- pred_21_all$ac_pred
lbc_target_21$ac_pred_filt <- pred_21_filt$ac_pred
lbc_target_21$ac_pred_all_filt <- pred_21_all_filt$ac_pred
lbc_target_21$ac_pred_dan <- pred_21_dan$ac_pred


## Correlation
####################################################################

# Non-log
r_36 <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred, use="pairwise.complete.obs") # 0.4444912
r_36_all <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_all, use="pairwise.complete.obs") # 0.405975
r_36_filt <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_filt, use="pairwise.complete.obs") # 0.4638506 (0.4446818 Carreras)
r_36_all_filt <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_all_filt, use="pairwise.complete.obs") # 0.470762 (0.4702639 Carreras)
r_36_dan <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_dan, use="pairwise.complete.obs") # 0.3624084

r_21 <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred, use="pairwise.complete.obs") # 0.4230204
r_21_all <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_all, use="pairwise.complete.obs") # 0.3802862
r_21_filt <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_filt, use="pairwise.complete.obs") # 0.4592133 (0.4359467 Carreras)
r_21_all_filt <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_all_filt, use="pairwise.complete.obs") # 0.4496676 (0.4422971 Carreras)
r_21_dan <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_dan, use="pairwise.complete.obs") # 0.2830282

# Log
r_36 <- cor(lbc_target_36$alcunitsupw_log, lbc_target_36$ac_pred, use="pairwise.complete.obs") # 0.3854889
r_36_all <- cor(lbc_target_36$alcunitsupw_log, lbc_target_36$ac_pred_all, use="pairwise.complete.obs") # 0.3480113
r_36_filt <- cor(lbc_target_36$alcunitsupw_log, lbc_target_36$ac_pred_filt, use="pairwise.complete.obs") # 0.4040634
r_36_all_filt <- cor(lbc_target_36$alcunitsupw_log, lbc_target_36$ac_pred_all_filt, use="pairwise.complete.obs") # 0.4151113
r_36_dan <- cor(lbc_target_36$alcunitsupw_log, lbc_target_36$ac_pred_dan, use="pairwise.complete.obs") # 0.3353729

r_21 <- cor(lbc_target_21$alcunitsupw_log, lbc_target_21$ac_pred, use="pairwise.complete.obs") # 0.3930249
r_21_all <- cor(lbc_target_21$alcunitsupw_log, lbc_target_21$ac_pred_all, use="pairwise.complete.obs") # 0.3552321
r_21_filt <- cor(lbc_target_21$alcunitsupw_log, lbc_target_21$ac_pred_filt, use="pairwise.complete.obs") # 0.4212163
r_21_all_filt <- cor(lbc_target_21$alcunitsupw_log, lbc_target_21$ac_pred_all_filt, use="pairwise.complete.obs") # 0.4089689
r_21_dan <- cor(lbc_target_21$alcunitsupw_log, lbc_target_21$ac_pred_dan, use="pairwise.complete.obs") # 0.2591051


## R2
####################################################################

# Incremental R2?
null_36 <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36 <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred, data=lbc_target_36))$r.squared
round(100*(full_36 - null_36), 3) # 19.259
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_all <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_all <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_36))$r.squared
round(100*(full_36_all - null_36_all), 3) # 16.137
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_36))$r.squared
round(100*(full_36_filt - null_36_filt), 3) # 20.838
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_all_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_all_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_36))$r.squared
round(100*(full_36_all_filt - null_36_all_filt), 3) # 21.39
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_dan <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_dan <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_36))$r.squared
round(100*(full_36_dan - null_36_dan), 3) # 12.475
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_36)) # p-value: < 2.2e-16

null_21 <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21 <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred, data=lbc_target_21))$r.squared
round(100*(full_21 - null_21), 3) # 18.311
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_all <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_all <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_21))$r.squared
round(100*(full_21_all - null_21_all), 3) # 14.167
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_21))$r.squared
round(100*(full_21_filt - null_21_filt), 3) # 22.311
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_all_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_all_filt <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_21))$r.squared
round(100*(full_21_all_filt - null_21_all_filt), 3) # 21.297
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_dan <- summary(lm(alcunitsupw ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_dan <- summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_21))$r.squared
round(100*(full_21_dan - null_21_dan), 3) # 7.544
summary(lm(alcunitsupw ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_21)) # p-value: 2.952e-15


# Incremental R2 (Log)?
null_36 <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36 <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred, data=lbc_target_36))$r.squared
round(100*(full_36 - null_36), 3) # 14.677
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_all <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_all <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_36))$r.squared
round(100*(full_36_all - null_36_all), 3) # 12.058
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_36))$r.squared
round(100*(full_36_filt - null_36_filt), 3) # 15.856
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_all_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_all_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_36))$r.squared
round(100*(full_36_all_filt - null_36_all_filt), 3) # 16.636
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_36)) # p-value: < 2.2e-16

null_36_dan <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_36))$r.squared
full_36_dan <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_36))$r.squared
round(100*(full_36_dan - null_36_dan), 3) # 10.641
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_36)) # p-value: < 2.2e-16

null_21 <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21 <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred, data=lbc_target_21))$r.squared
round(100*(full_21 - null_21), 3) # 16.031
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_all <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_all <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_21))$r.squared
round(100*(full_21_all - null_21_all), 3) # 12.49
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_filt, data=lbc_target_21))$r.squared
round(100*(full_21_filt - null_21_filt), 3) # 19.024
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_all_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_all_filt <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_21))$r.squared
round(100*(full_21_all_filt - null_21_all_filt), 3) # 17.872
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_all_filt, data=lbc_target_21)) # p-value: < 2.2e-16

null_21_dan <- summary(lm(alcunitsupw_log ~ age + as.factor(sex), data=lbc_target_21))$r.squared
full_21_dan <- summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_21))$r.squared
round(100*(full_21_dan - null_21_dan), 3) # 6.31
summary(lm(alcunitsupw_log ~ age + as.factor(sex) + ac_pred_dan, data=lbc_target_21)) # p-value: 2.952e-15

